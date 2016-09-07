/*
 * PhyloTree.h
 *	An binary Phylogenic tree class
 *	the tree is rooted if the root node has a null parent, or unrooted if the parent is actually the 3'rd child
 *  Created on: Mar 25, 2016
 *      Author: zhengqi
 */

#ifndef SRC_PHYLOTREE_H_
#define SRC_PHYLOTREE_H_

#include <string>
#include <vector>
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

#include "HmmUFOtuConst.h"
#include "StringUtils.h"
#include "SeqCommons.h"
#include "DigitalSeq.h"
#include "MSA.h"

namespace EGriceLab {

using std::string;
using std::vector;
using std::istream;
using std::ostream;
using Eigen::Matrix4Xd;

<<<<<<< HEAD
struct PhyloTree;
=======
class PhyloTree {
	/* nested types and enums */
private:
	class PhyloTreeNode;
	typedef PhyloTreeNode PTNode;

	class PhyloTreeNode {
	public:
		/* constructors */
		/* Default constructor */
		PhyloTreeNode() : parent(NULL), childL(NULL), childR(NULL),
		length(0) { }
>>>>>>> refs/heads/master

<<<<<<< HEAD
typedef PhyloTree PT;
=======
		/* construct a node with a given DigitalSeq */
		explicit PhyloTreeNode(const DigitalSeq& seq) :
				seq(seq), parent(NULL), childL(NULL), childR(NULL),
				length(0) { }
>>>>>>> refs/heads/master

<<<<<<< HEAD
//typedef boost::variant<boost::recursive_wrapper<PT> > pt_node;
=======
		/* construct a node with a given PrimarySeq */
		explicit PhyloTreeNode(const PrimarySeq& seq) :
				seq(seq), parent(NULL), childL(NULL), childR(NULL),
				length(0) { }
>>>>>>> refs/heads/master

<<<<<<< HEAD
struct PhyloTree {
=======
		virtual ~PhyloTreeNode() { }

		/* Member methods */
		bool hasParent() const {
			return parent != NULL;
		}

		bool hasChildL() const {
			return childL != NULL;
		}

		bool hasChildR() const {
			return childR != NULL;
		}

		PTNode* getParent() const {
			return parent;
		}

		const PTNode* getParent() const {
			return parent;
		}

		PTNode* getChildL() const {
			return childL;
		}

		const PTNode* getChildL() const {
			return childL;
		}

		PTNode* getChildR() const {
			return childR;
		}

		const PTNode* getChildR() const {
			return childR;
		}

		/** test whether this node is a root node */
		bool isRoot() const;

		/** test whether this node is a leaf node */
		bool isLeaf() const;

		/** test whether this node is a tip node */
		bool isTip() const;

		DigitalSeq seq;

		PhyloTreeNode* parent;
		PhyloTreeNode* childL;
		PhyloTreeNode* childR;

		double length;

		Matrix4Xd cost; /* cost (negative logLiklihood) of observing this sequence given the mode and the tree */
	};

public:
>>>>>>> refs/heads/master
	/* constructors */
	/* Default constructor */
	PhyloTree() : length(0) {
		/* Assert IEE559 at construction time */
		assert(std::numeric_limits<double>::is_iec559);
	}

	/* Member methods */

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

	/** test whether the cost is evaluated */
	bool isEvaluated() const {
		return cost.cols() > 0 && (cost.array() != inf).any();
	}

	/** Get the number of aligned sites of this tree */
	int alnSites() const;

	/** Get the number of nodes of this tree using Dfs search */
	int numNodes() const;

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

	/* friend operators */
	friend ostream& operator<<(ostream& out, const PT& tree);

	/* member fields */
	string name;
	double length; /* branch length of this node */
	vector<PT> children;

	DigitalSeq seq;
	Matrix4Xd cost; /* cost (negative log liklihood) of observing this sequence given the model and the tree */

};

<<<<<<< HEAD
inline int PhyloTree::readTree(const string& treefn, const string& format, const MSA* msa) {
	/* Read tree structure */
	int nNodes = readTree(treefn, format);
	if(nNodes == -1) {
		std::cerr << "Failed to read in tree file: " << treefn << " in " << format << " format" << std::endl;
		return nNodes;
=======
inline bool PhyloTree::PhyloTreeNode::isRoot() const {
	return (hasParent() && hasChildL() && hasChildR()) /* internally rooted */ ||
			(isLeaf() && hasParent()) /* leaf rooted */;
}

inline bool PhyloTree::PhyloTreeNode::isLeaf() const {
	return childL == NULL && childR == NULL;
}

inline bool PhyloTree::PhyloTreeNode::isTip() const {
	return !isLeaf() && childL->isLeaf() && childR->isLeaf();
}

inline void PhyloTree::swap(PhyloTree& other) {
	using std::swap;
	swap(root, other.root);
	swap(abc, other.abc);
}

inline bool PhyloTree::isRooted() const {
	return root->isRoot();
}

inline int PhyloTree::numSites() const {
	return root->seq.length();
}

inline int PhyloTree::numNodes() const {
	return dfsNodes().size();
}

inline bool PhyloTree::isSibling(const PhyloTreeNode* node1, const PhyloTreeNode* node2) const {
	assert(isRooted()); // only rooted tree can test siblings
	return !node1->isRoot() && !node2->isRoot() && node1->parent == node2->parent;
}

inline set<const PhyloTree::PhyloTreeNode*> PhyloTree::children(const PhyloTree::PhyloTreeNode* node) const {
	set<const PhyloTree::PhyloTreeNode*> children;
	children.insert(node->childL);
	children.insert(node->childR);
	return children;
}

inline set<const PhyloTree::PhyloTreeNode*> PhyloTree::ancestors(const PhyloTree::PhyloTreeNode* node) const {
	set<const PhyloTree::PhyloTreeNode*> ancestors;
	while(!node->isRoot()) {
		ancestors.insert(node->parent);
		node = node->parent;
>>>>>>> refs/heads/master
	}
	/* Read MSA sequences */
	int nSeqs = readSeqFromMSA(msa);
	if(nSeqs == -1) {
		std::cerr << "Failed to read in the MSA" << std::endl;
		return nSeqs;
	}
	return nNodes;
}

inline int PhyloTree::readTree(const string& treefn, const string& format) {
	if(StringUtils::toLower(format) == "newick")
		return readNiwickTree(treefn);
	else
		throw invalid_argument("Unsupported tree file format '" + format + "'");
}

} /* namespace EGriceLab */

#endif /* SRC_PHYLOTREE_H_ */

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
	newick_grammar() : newick_grammar::base_type(tree) {
		using phoenix::at_c;
		using phoenix::push_back;

		/* label grammars */
		/* unquoted label is printable characters
		 * without blanks, parentheses, brackets, single_quotes, colons, semicolons or commas
		 */
		unquoted_label %= qi::lexeme[+(ascii::print - ascii::space - '(' - ')' - '[' - ']' - '\'' - ':' - ';' - ',') ];
		/* quoted label is ' printable characters ' */
		quoted_label %= '\'' >> qi::lexeme[+(ascii::print)] >> '\'';

		label = unquoted_label | quoted_label;

		/* branch length grammar */
		branch_length %= ':' >> qi::double_;

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
	}

	private:
	qi::rule<Iterator, PT()> tree;
	qi::rule<Iterator, PT()> subtree;
	qi::rule<Iterator, std::vector<PT>()> descendant_list;
	qi::rule<Iterator, double()> branch_length;
	qi::rule<Iterator, std::string()> unquoted_label, quoted_label, label;
};

inline ostream& operator<<(ostream& out, const PT& tree) {
	bool first = true;
	if(!tree.children.empty()) {
		out << '(';
		for(std::vector<PT>::const_iterator it = tree.children.begin(); it != tree.children.end(); ++it) {
			out << (first ? "" : ",") << *it;
			first = false;
		}
		out << ')';
	}
	if(!tree.name.empty())
		out << tree.name;
	if(tree.length > 0)
		out << ':' << tree.length;

	out << ';';
	return out;
}

inline int PhyloTree::alnSites() const {
	return seq.length();
}

} /* namespace EGriceLab */
