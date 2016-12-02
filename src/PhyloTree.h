/*
 * PhyloTree.h
 *	An binary Phylogenic tree class
 *	the tree is rooted if the root node has a null parent, or unrooted if the parent is actually the 3'rd child
 *  Created on: Mar 25, 2016
 *      Author: zhengqi
 */

#ifndef PHYLOTREE_H_
#define PHYLOTREE_H_

#include <string>
#include <vector>
#include <set>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <cstddef>
#include <cassert>
#include <eigen3/Eigen/Dense>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/variant/recursive_variant.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include "HmmUFOtuConst.h"
#include "StringUtils.h"
#include "SeqCommons.h"
#include "DigitalSeq.h"
#include "MSA.h"

namespace EGriceLab {

using std::string;
using std::vector;
using std::set;
using std::istream;
using std::ostream;
using Eigen::Matrix4Xd;
using Eigen::RowVectorXd;

struct PhyloTree;

typedef PhyloTree PT;

//typedef boost::variant<boost::recursive_wrapper<PT> > pt_node;

struct PhyloTree {
	/** Nested types and enums */
//	enum TaxonRank { K, P, C, O, F, G, S };

	/*
	 * Pre-calculated PhyloTree new node placement cost at given node's branch, considering all possible alphabet symbols
	 */
	struct NodePlacementCost {
		/* constructors */
		/** default constructor */
		NodePlacementCost() : id(-1), annoDist(0), K(0), L(0)
		{ 	}

		/* member methods */
		istream& load(istream& in);
		ostream& save(ostream& out) const;

		/* member fields */
		long id;
		string anno;
		double annoDist;
		int K; /* number of alphabet size + 1 */
		int L; /* number of aligned sites */
		MatrixXd treeCost; /* placement cost matrix with K x L elements indicating pre-calculated final cost of the entire tree */
	};

	/* constructors */
	/* Default constructor */
	PhyloTree() : length(0), id(-1), annoDist(0), parent(NULL) {
		/* Assert IEE559 at construction time */
		assert(std::numeric_limits<double>::is_iec559);
	}

	/* Member methods */

	/** test whether this node is named */
	bool isNamed() const {
		return !name.empty();
	}

	/** test whether this SubTree is a leave root */
	bool isLeafRoot() const {
		return children.size() == 1;
	}

	/** test whether this SubTree is an internal root */
	bool isInternalRoot() const {
		return children.size() > 2;
	}

	/** test whether this SubTree is a root node */
	bool isRoot() const {
		return isLeafRoot() || isInternalRoot();
	}

	/** test whether this node is a leaf node */
	bool isLeaf() const {
		return children.empty();
	}

	/** test whether this node is a tip node */
	bool isTip() const {
		if(isLeaf())
			return false;
		/* test whether all children are leaves */
		for(vector<PT>::const_iterator it = children.begin(); it != children.end(); ++it)
			if(!it->isLeaf())
				return false;
		return true;
	}

	/** test whether this tree node is evaluated */
	bool isEvaluated() const {
		return cost.cols() > 0 && (cost.array() != inf).any();
	}

	/** test whether this tree is evaluated at specific site */
	bool isEvaluated(int j) const {
		return cost.cols() > 0 && (cost.col(j).array() != inf).any();
	}

	/** reset the (evaluated) cost of this tree */
	void resetCost() {
		if(cost.cols() > 0)
			cost.setConstant(inf);
		if(scale.rows() > 0)
			scale.setZero();
	}

	/**
	 * Get the scaled cost of this node
	 */
	Matrix4Xd getCost() const {
		return cost.rowwise() + scale;
	}

	/**
	 * Re-scale the cost at position j of this node by subtracting a constant
	 */
	void reScale(double c, int j) {
		cost.col(j).array() -= c;
		scale(j) += c;
	}

	/**
	 * Re-scale the cost of all position by subtracting a constant
	 */
	void reScale(double c) {
		cost.array() -= c;
		scale.array() += c;
	}

	/**
	 * Re-scale the cost of all position by subtracting a given vector
	 */
	void reScale(const RowVectorXd& c) {
		cost.rowwise() -= c;
		scale += c;
	}

	/** Get the number of aligned sites of this tree */
	int alnSites() const;

	/** Get the number of nodes of this tree using Dfs search */
	int numNodes() const;

	/** Get the number of leaf node of this tree using Dfs search */
	int numLeaves() const;

	/**
	 * Get a set of const pointers of all leaf nodes of this sub-tree
	 */
	set<const PT*> leafNodes() const;

	/**
	 * Get the first leave node of this sub-tree
	 */
	const PT* firstLeaf() const;

	/**
	 * Read the tree structure and sequence from an input file of given format
	 * @param treefn  tree filename
	 * @param format  tree file format
	 * @param msa  Multiple Sequence Alignment of this tree
	 * @throw illegal_argument exception if is not a supported file format
	 */
	int readTree(const string& treefn, const string& format, const MSA* msa);

	/**
	 * Read the tree structure (but not the sequence) from an input file of given format
	 * @param treefn  tree filename
	 * @param format  tree file format
	 * @param msa  Multiple Sequence Alignment of this tree
	 * @return  number of tree nodes read, or -1 if anything failed
	 * @throw illegal_argument exception if is not a supported file format
	 */
	int readTree(const string& treefn, const string& format);

	/**
	 * Read the tree structure (but not the sequence) from a newick file
	 * @param treefn  tree filename
	 * @return
	 */
	int readNiwickTree(const string& treefn);

	/**
	 * Read the sequence in the this tree according to given MSA
	 * @return  sequences with MSA assigned, or -1 if anything failed
	 */
	int readSeqFromMSA(const MSA* msa);

	/**
	 * print this tree in newick format
	 * @param out  output stream
	 * return  output stream modified
	 */
	ostream& print(ostream& out) const;

	/**
	 * Set the unique ID and parent pointer for this sub-tree recursively using DFS algorithm,
	 * @param initID  starting ID for the root node
	 * @return  total nodes found in this sub-tree
	 */
	long setIDandParent();

	/**
	 * update the parent pointer for this sub-tree recursively using DFS algorithm
	 * so nodes will not point to parents outside the current tree. Should be called whenever the tree is copied
	 */
	long updateParent();

	/**
	 * Annotate this sub-tree recursively using DFS algorithm,
	 * where unnamed nodes will be annotated to their closest named ancestor node
	 * and named node will have their annotation as their name
	 */
	void annotate();

	/* friend operators */
	friend ostream& operator<<(ostream& out, const PT& tree);

	/* member fields */
	long id; /* unique id for each node, in Depth-first search (DFS) order */
	string name; /* node name */
	double length; /* branch length of this node */
	string anno; /* node annotation */
	double annoDist; /* Phylogenetic distance to the named ancestor node from which this node get its annotation */
	vector<PT> children;
	const PT* parent; /* convenient pointer to this PTNode's parent */

	DigitalSeq seq;
	Matrix4Xd cost; /* cost (negative log likelihood) of observing this sequence given the model and the tree */
	RowVectorXd scale; /* log-scale constant used during calculating and storing the cost to avoid numeric underfly */

	/* static fields */
	static const double MIN_EXPONENT;
};

inline int PhyloTree::readTree(const string& treefn, const string& format, const MSA* msa) {
	/* Read tree structure */
	int nNodes = readTree(treefn, format);
	if(nNodes == -1) {
		std::cerr << "Failed to read in tree file: " << treefn << " in " << format << " format" << std::endl;
		return nNodes;
	}
	int nLeaves = numLeaves();
	std::cerr << "PhyloTree read with " << nNodes << " nodes and " << nLeaves << " leaves" << std::endl;
	/* Read MSA sequences */
	int nSeqs = readSeqFromMSA(msa);
	std::cerr << "Read in " << nSeqs << " sequences from MSA" << std::endl;
	if(nSeqs == -1) {
		std::cerr << "Failed to read in the MSA" << std::endl;
		return nSeqs;
	}
	else if(nSeqs != numLeaves()) {
		std::cerr << "Unmatched leaves in tree file " << treefn << " and MSA file " << msa << std::endl;
		return -1;
	}
	else
		return nNodes;
}

inline int PhyloTree::readTree(const string& treefn, const string& format) {
	if(StringUtils::toLower(format) == "newick")
		return readNiwickTree(treefn);
	else {
		std::cerr << "Unsupported tree file format '" << format << "'" << std::endl;
		return -1;
	}
}

} /* namespace EGriceLab */


// adapt the structure to fusion phoenix
BOOST_FUSION_ADAPT_STRUCT(
	EGriceLab::PhyloTree,
	(std::string, name)
	(double, length)
	(std::vector<EGriceLab::PT>, children)
	(EGriceLab::PT*, parent)
	(const EGriceLab::DigitalSeq*, seq)
	(Eigen::Matrix4Xd, cost)
)

namespace EGriceLab {
/* namespace aliasing */
namespace qi = boost::spirit::qi;
namespace phoenix = boost::phoenix;
namespace fusion = boost::fusion;
namespace ascii = boost::spirit::ascii;

/* This generic grammar parse a generic iterator to a PhyloTree */
template <typename Iterator>
struct newick_grammar :
		qi::grammar<Iterator, PT()> {
	/* grammar constructor */
	newick_grammar() : newick_grammar::base_type(tree, "Newick tree") {
		using ascii::char_;
		using ascii::string;
		using phoenix::at_c;
		using phoenix::push_back;
		using qi::on_error;
		using qi::fail;

		using phoenix::construct;
		using phoenix::val;

		/* label grammars */
		/* unquoted label is printable characters
		 * without blanks, parentheses, brackets, single_quotes, colons, semicolons or commas
		 */
		unquoted_label %= qi::lexeme[+(ascii::print - ascii::space - '(' - ')' - '[' - ']' - '\'' - ':' - ';' - ',') ];
		/* quoted label is ' printable characters ' */
		quoted_label %= '\'' > qi::lexeme[+(ascii::print - '\'')] > '\'';

		label = unquoted_label | quoted_label;

		/* branch length grammar */
		branch_length %= ':' > qi::double_;

		/* subtree grammar */
		subtree =
			// assign vector of children to the third element of PT, optional
			-(descendant_list  [at_c<2>(qi::_val) = qi::_1])
			// assign the label to the first element, optional
			>> -(label  [at_c<0>(qi::_val) = qi::_1])
			// assign the branch length to the second element, optional
			>> -(branch_length  [at_c<1>(qi::_val) = qi::_1]);

		/* descentdant_list is a vector of PT, which will be auto pushed to the compatible PT vector */
		descendant_list %=
			'(' >> subtree % ',' >> ')'; /* this grammar has the attribute of vector<PT> */

		// The tree receive the whole subtree using %=
		tree %= subtree >> ';';

		unquoted_label.name("unquoted_label");
		quoted_label.name("quoted_label");
		label.name("label");
		branch_length.name("branch_length");
		subtree.name("subtree");
		descendant_list.name("descendant_list");
		tree.name("tree");

		on_error<fail>(
				tree,
				std::cout
					<< val("Error! Expecting ")
					<< qi::_4                             // what failed?
					<< val(" here:\"")
					<< construct<std::string>(qi::_3, qi::_2) // iterators to error-pos, end
					<< val("\"")
					<< std::endl
				);
	}

	private:
	qi::rule<Iterator, PT()> tree;
	qi::rule<Iterator, PT()> subtree;
	qi::rule<Iterator, std::vector<PT>()> descendant_list;
	qi::rule<Iterator, double()> branch_length;
	qi::rule<Iterator, std::string()> unquoted_label, quoted_label, label;
};

inline ostream& operator<<(ostream& out, const PT& tree) {
	tree.print(out);
	return out << ';';
}

inline int PhyloTree::alnSites() const {
	return seq.length();
}

} /* namespace EGriceLab */

#endif /* PHYLOTREE_H_ */
