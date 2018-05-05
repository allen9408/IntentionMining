package ICHI2018;



/*
 * This method is based on the maximum likelihood
 * User can predetermine the state number or just let the program run to a point 
 * that furthur splitting won't increaset the probability any more
 * 
 * Written by Sen Yang 

 */


import java.util.ArrayList;
import java.util.List;
import ICHI2018.HierarchicalHmm;

public class HierarchicalHmmList {
	private int depth = 0;
	private int size = 0;
	private List<HierarchicalHmm> HmmList;

	public HierarchicalHmmList(int depth) {
		this.HmmList = new ArrayList<HierarchicalHmm>();
		this.depth = depth;
	}

	public void addHmm(HierarchicalHmm hmm) {
		this.HmmList.add(hmm);
		this.size++;
	}

	public List<HierarchicalHmm> getHmmList() {
		return this.HmmList;
	}

	public int getDepth() {
		return this.depth;
	}
	public int size() {
		return this.size;
	}
	public HierarchicalHmm getHmmByIndex(int index) {
		return this.HmmList.get(index);
	}
	public List<HierarchicalHmm> getHmmByParent(int parent) {
		List<HierarchicalHmm> result = new ArrayList<HierarchicalHmm>();
		for (HierarchicalHmm hmm : this.HmmList) {
			if (hmm.getParent() == parent) {
				result.add(hmm);
			}
		}
		return result;
	}
	public int indexOf(HierarchicalHmm hmm) {
		int index = -1;
		for (int i = 0; i < this.size; i++) {
			HierarchicalHmm tmpHmm = this.HmmList.get(i);
			if ((hmm.getDepth() == tmpHmm.getDepth()) && (hmm.getParent() == tmpHmm.getParent()) && (hmm.getStateNum() == tmpHmm.getStateNum()))
				return i;
		}
		return index;
	}
}
