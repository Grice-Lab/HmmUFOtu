/*
 * PhyloTreeUnrooted.h
 *  An Unrooted Phylogenic Tree (PTUnrooted)
 *  A PTUnrooted can be rooted at any node
 *  A strict PTUnrooted node has exactly three (3) neighbors of its internal nodes
 *  When rooted, at most 1 neighbor is designated as its parent
 *  Tree nodes are number indexed internally
 *  Created on: Dec 1, 2016
 *      Author: zhengqi
 */

#ifndef SRC_PHYLOTREEUNROOTED_H_
#define SRC_PHYLOTREEUNROOTED_H_

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
#include <boost/shared_ptr.hpp>

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
using boost::shared_ptr;

class PhyloTreeUnrooted; /* forward declaration */

typedef PhyloTreeUnrooted PTUnrooted;

class PhyloTreeUnrooted {
public:
	/* nested types and enums */

	struct PhyloTreeUnrootedNode;
	typedef PTUnrooted::PhyloTreeUnrootedNode PTUNode;

	typedef shared_ptr<PTUNode> PTUNodePtr; /* use boost shared_ptr to hold node pointers */
	typedef shared_ptr<const PTUNode> PTUNodeConstPtr; /* a const version of boost shared_ptr */

	/**
	 * A PTUnrooed node that stores its neighbor pointers
	 * if the tree is rooted, it always position 0 to store the parent
	 * and position 1..N to store the children
	 */
	struct PhyloTreeUnrootedNode {
		/* constructors */
		/**
		 * Default constructor, constructing a node with parent set to NULL
		 * and parent length set to 0
		 * and parent cost not initiated
		 */
		PhyloTreeUnrootedNode() : hasParent(false) {
			neighbors.push_back(PTUNodePtr());
			length.push_back(0.0);
			cost.push_back(Matrix4Xd());
		}

		/* Member methods */
		/** test whether this node is initiated */
		bool isInit() const {
			return !neighbors.empty();
		}

		/** test whether this node is named */
		bool isNamed() const {
			return !name.empty();
		}

		/** test whether this is a leave node */
		bool isLeaf() const {
			return neighbors.size() == 1;
		}

		/** test whether this is an internal node */
		bool isInternal() const {
			return neighbors.size() > 1;
		}

		/** test whether this is a leaf root */
		bool isLeafRoot() const {
			return isLeaf() && !hasParent;
		}

		/** test whether this is an internal root */
		bool isInternalRoot() const {
			return isInternal() && !hasParent;
		}

		/** test whether this is a root */
		bool isRoot() const {
			return !hasParent;
		}

		/** test whether this node is a tip node */
		bool isTip() const {
			if(neighbors.empty())
				return false;
			/* test whether all children are leaves, ignore the parent */
			for(vector<PTUNodePtr>::const_iterator child = neighbors.begin() + 1; child != neighbors.end(); ++child)
				if(!(*child)->isLeaf())
					return false;
			return true;
		}

		/** reset the incoming cost of this node */
		void resetCost();

		/** reset the incoming cost of a neighbor */
		void resetCost(int i) {
			if(cost[i].cols() > 0)
				cost[i].setConstant(inf);
		}

		/**
		 * Get the cost of this node be merging incoming costs
		 */
		Matrix4Xd getCost() const;

		string name; /* node name, need to be unique for database loading */
		DigitalSeq seq; /* sequence of this node */
		vector<PTUNodePtr> neighbors; /* pointers to neighbors */
		vector<double> length; /* branch lengths to neighbors */
		bool hasParent; /* flag used when this node has a parent, which is always the first neighbor */
		vector<Matrix4Xd> cost; /* incoming cost (message) from its children roots, when the tree is rooted */
	};

	PhyloTreeUnrooted() {
		// TODO Auto-generated constructor stub

	}

	/** Get the number of nodes of this tree using Dfs search */
	int numNodes() const;

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

private:
	int L; /* number of aligned sites */

	PTUNodePtr root; /* root node of this tree */
	vector<PTUNodePtr> nodes; /* indexed tree nodes */
	vector<string> anno; /* index node phylogenetic annotations */
	vector<double> annoDist; /* index node annotation distance */

public:
	/* static fields */
	static const int INTERNAL_NEIGHBORS = 3; /* a PTUnrooted internal node always has 3 neighbors */
	static const double MIN_EXPONENT;
};

inline void PhyloTreeUnrooted::PhyloTreeUnrootedNode::resetCost() {
	assert(neighbors.size() == cost.size());
	for(vector<PTUNodePtr>::size_type i = 0; i != neighbors.size(); ++i)
		resetCost(i);
}

inline int PhyloTreeUnrooted::readTree(const string& treefn, const string& format) {
	if(StringUtils::toLower(format) == "newick")
		return readNiwickTree(treefn);
	else {
		std::cerr << "Unsupported tree file format '" << format << "'" << std::endl;
		return -1;
	}
}


} /* namespace EGriceLab */

//// adapt the PTUNode to fusion phoenix
//BOOST_FUSION_ADAPT_STRUCT(
//	EGriceLab::PhyloTreeUnrooted::PhyloTreeUnrootedNode,
//	(std::string, name)
//	(EGriceLab::DigitalSeq, seq)
//	(std::vector<EGriceLab::PhyloTreeUnrooted::PTUNodePtr>, neighbors)
//	(std::vector<double>, length)
//	(bool, hasParent)
//	(std::vector<Eigen::Matrix4Xd>, cost)
//)

// contruct a PTUNode parser
namespace EGriceLab {
/* namespace aliasing */
namespace qi = boost::spirit::qi;
namespace phoenix = boost::phoenix;
namespace fusion = boost::fusion;
namespace ascii = boost::spirit::ascii;

/* This generic grammar parse a generic iterator to a PTUNodePtr */
template <typename Iterator>
struct newick_grammar :
		qi::grammar<Iterator, PTUnrooted::PTUNodePtr()> {
	/* grammar constructor */
	newick_grammar() : newick_grammar::base_type(tree, "Newick tree") {
		using ascii::char_;
		using ascii::string;
		using phoenix::at_c;
		using phoenix::push_back;
		using qi::eps;
		using qi::_val;
		using qi::on_error;
		using qi::fail;

		using phoenix::construct;
		using phoenix::ref;
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
		subtree = eps [qi::_val = phoenix::new_<PTUnrooted::PTUNode> ] >> /* allocate the node */
			// assign vector of children to the third element of PT, optional
			-(descendant_list  ) /* get the neighbors, optional */
			// assign the label, optional
			>> -(label  )
			// assign the branch length, optional
			>> -(branch_length  );

		/* descentdant_list is a vector of PTUNodePtr, which will be auto pushed to the compatible PTUNodePtr vector */
		descendant_list %=
			'(' >> subtree % ',' >> ')'; /* this grammar has the attribute of vector<PTUNodePtr> */

		// The tree receive the whole subtree using %=
		tree %= subtree >> ';';

		unquoted_label.name("unquoted_label");
		quoted_label.name("quoted_label");
		label.name("label");
		branch_length.name("branch_length");
		subtree.name("subtree");
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
	qi::rule<Iterator, PTUnrooted::PTUNodePtr()> tree;
	qi::rule<Iterator, PTUnrooted::PTUNodePtr()> subtree;
	qi::rule<Iterator, std::vector<PTUnrooted::PTUNodePtr>()> descendant_list;
	qi::rule<Iterator, double()> branch_length;
	qi::rule<Iterator, std::string()> unquoted_label, quoted_label, label;

};

} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREEUNROOTED_H_ */
