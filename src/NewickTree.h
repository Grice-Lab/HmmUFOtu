/*******************************************************************************
 * This file is part of HmmUFOtu, an HMM and Phylogenetic placement
 * based tool for Ultra-fast taxonomy assignment and OTU organization
 * of microbiome sequencing data with species level accuracy.
 * Copyright (C) 2017  Qi Zheng
 *
 * HmmUFOtu is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HmmUFOtu is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AlignerBoost.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/*
 * NewickTree.h
 *  A NewickTree class
 *  A NewickTree is a rooted, bi/multi-furcating phylogenetic tree
 *  A NewickTree only stores the basic relationship between tree nodes and their children, and branch length to parent
 *  but not their parent directly
 *  Tree nodes can be unnamed or even duplicated
 *  Created on: Dec 2, 2016
 *      Author: zhengqi
 */

#ifndef SRC_NEWICKTREE_H_
#define SRC_NEWICKTREE_H_

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

namespace EGriceLab {

using std::string;
using std::vector;
using std::istream;
using std::ostream;

struct NewickTree;

typedef NewickTree NT;

struct NewickTree {
public:
	/* constructors */
	/** Default constructor */
	NewickTree() : length(0) { }

	/** Construct a Newick tree node with given name and an optional parent distance */
	explicit NewickTree(const string& name, double length = 0) : name(name), length(length)
	{ }

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

	/** test whether this node is a leaf node */
	bool isLeaf() const {
		return children.empty();
	}

	/** remove all offspring nodes of this subtree */
	void clear() {
		children.clear();
	}

	/** add a child to this NT node */
	void addChild(const NT& node) {
		children.push_back(node);
	}

	/**
	 * Read the tree structure from input in Newick format
	 * @param in  input stream
	 * @return  the modified input
	 */
	istream& read(istream& in);

	/**
	 * Write the tree structure to a output in Newick format
	 * @param out  output stream
	 * @return  the modified output
	 */
	ostream& write(ostream& out) const;

	/* non-member functions */
	friend istream& operator>>(istream& in, NT& tree);

	friend ostream& operator<<(ostream& out, const NT& tree);

	/* member fields */
	string name; /* subtree (node) name */
	double length; /* branch length (to parent) of this subtree */
	vector<NT> children;

	/* static fields */
	static const string& INVALID_CHARS;
}; /* struct NewickTree */

inline istream& operator>>(istream& in, NT& tree) {
	tree.read(in);
	return in;
}

inline ostream& operator<<(ostream& out, const NT& tree) {
	tree.write(out);
	return out << ';';
}

} /* namespace EGriceLab */

// adapt the structure to fusion phoenix
BOOST_FUSION_ADAPT_STRUCT(
	EGriceLab::NewickTree,
	(std::string, name)
	(double, length)
	(std::vector<EGriceLab::NT>, children)
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
		qi::grammar<Iterator, NT()> {
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
			// assign vector of children to the third element of NT, optional
			-(descendant_list  [at_c<2>(qi::_val) = qi::_1])
			// assign the label to the first element, optional
			>> -(label  [at_c<0>(qi::_val) = qi::_1])
			// assign the branch length to the second element, optional
			>> -(branch_length  [at_c<1>(qi::_val) = qi::_1]);

		/* descentdant_list is a vector of NT, which will be auto pushed to the compatible NT vector */
		descendant_list %=
			'(' >> subtree % ',' >> ')'; /* this grammar has the attribute of vector<NT> */

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
	qi::rule<Iterator, NT()> tree; /* tree is a NT node */
	qi::rule<Iterator, NT()> subtree; /* so is subtree */
	qi::rule<Iterator, std::vector<NT>()> descendant_list;
	qi::rule<Iterator, double()> branch_length;
	qi::rule<Iterator, std::string()> unquoted_label, quoted_label, label;
};

} /* namespace EGriceLab */

#endif /* SRC_NEWICKTREE_H_ */
