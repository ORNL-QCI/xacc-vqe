
#ifndef VQE_TRANSFORMATION_BK_FENWICK_HPP_
#define VQE_TRANSFORMATION_BK_FENWICK_HPP_

#include <vector>

class FNode : public std::enable_shared_from_this<FNode> {

public:
	using FNodePtr = std::shared_ptr<FNode>;

	FNodePtr parent;
	std::set<FNodePtr> children;
	int index;

	FNode(FNodePtr p, std::set<FNodePtr> cs) : index(-1), parent(p), children(cs) {}
	FNode(std::set<FNodePtr> cs) : index(-1), children(cs) {}

	std::set<FNodePtr> getAncestors() {

		auto node = shared_from_this();
		std::set<FNodePtr> ancestors;
		while (node->parent) {
			ancestors.insert(node);
			node = node->parent;
		}

		return ancestors;
	}
};

class FTree {

	using FNodePtr = std::shared_ptr<FNode>;

protected:

	std::vector<FNodePtr> nodes;

	FNodePtr root;
public:

	FTree(const int nQubits) {

		for (int i = 0; i < nQubits; i++) {
			nodes.push_back(std::make_shared<FNode>(std::set<FNodePtr>{}));
		}

		if (nQubits > 0 ) {
			root = nodes[nQubits-1];
			root->index = nQubits - 1;

			std::function<void(const int, const int, FNodePtr)> construct;
			construct =
					[this,&construct](const int idx1, const int idx2, FNodePtr parent) {
						if( idx1 >= idx2) {
							return;
						} else {
							auto pivot = (idx1 + idx2) >> 1;
							auto child = nodes[pivot];

							child->index = pivot;
							parent->children.insert(child);
							child->parent = parent;

							construct(idx1, pivot, child);
							construct(pivot + 1, idx2, parent);
						}
					};
			construct(0, nQubits-1, root);
		}
	}

	std::set<FNodePtr> getUpdateSet(const int i) {
		return nodes[i]->getAncestors();
	}

	std::set<FNodePtr> getChildrenSet(const int i) {
		return nodes[i]->children;
	}

	std::set<FNodePtr> getRemainderSet(const int i) {
		std::set<FNodePtr> rset;
		auto ancestors = getUpdateSet(i);

		for (auto a : ancestors) {
			for (auto c : a->children) {
				if (c->index < i) {
					rset.insert(c);
				}
			}
		}

		return rset;
	}

	std::set<FNodePtr> getParitySet(const int i) {
		auto result = getRemainderSet(i);
		auto cs = getChildrenSet(i);

		std::set<FNodePtr> paritySet = result;
		paritySet.insert(cs.begin(), cs.end());

		return paritySet;
	}
};


#endif /* VQE_TRANSFORMATION_BK_FENWICK_HPP_ */
