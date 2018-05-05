package ICHI2018;

import java.util.ArrayList;
import java.util.List;
import ICHI2018.HierarchicalHmm;
import ICHI2018.HierarchicalHmmList;
import HMM.DataType.Hmm;
import HMM.Observation.ObservationDiscrete;
import ICHI2018.SampleData;

public class HierarchicalHmmSet {
	private List<HierarchicalHmmList> hmmSet;

	public HierarchicalHmmSet() {
		this.hmmSet = new ArrayList<HierarchicalHmmList>();
	}

	public int depth() {
		return this.hmmSet.size();
	}

	public void clear() {
		this.hmmSet.clear();
	}
	public int size() {
		return this.hmmSet.size();
	}

	public int levelSize(int level) {
		if (level > this.hmmSet.size() - 1) {
			System.out.println("Get level: " + level + ", out of: " + this.hmmSet.size());
			return 0;
		}
		return this.hmmSet.get(level).size();

	}

	public List<Hmm<ObservationDiscrete<SampleData>>> getLevelList(int level) {
		if (level > this.hmmSet.size() - 1) {
			System.out.println("Get level: " + level + ", out of: " + this.hmmSet.size());
			return null;
		}
		List<HierarchicalHmm> hHmmList = hmmSet.get(level).getHmmList();
		List<Hmm<ObservationDiscrete<SampleData>>> hmmList = new ArrayList<Hmm<ObservationDiscrete<SampleData>>>();
		for (HierarchicalHmm hHmm : hHmmList) {
			hmmList.add(hHmm.getHmm());
		}
		return hmmList;
	}

	public void addHmm(Hmm<ObservationDiscrete<SampleData>> hmm, int depth, int stateNum, int parent) {
		HierarchicalHmm hHmm = new HierarchicalHmm(hmm, depth, stateNum, parent);
		// find if the depth exists
		int index = 0;
		while (depth > hmmSet.size()-1) {
			HierarchicalHmmList newList = new HierarchicalHmmList(index);
			this.hmmSet.add(newList);
			index++;
		}
		HierarchicalHmmList curLevel = hmmSet.get(depth);
		curLevel.addHmm(hHmm);
	}

	public Hmm<ObservationDiscrete<SampleData>> getHmmFromSet(int depth, int stateNum, int parent) {
		List<HierarchicalHmm> curLevel = this.hmmSet.get(depth).getHmmList();
		for (HierarchicalHmm hHmm : curLevel) {
			if (hHmm.getStateNum() == stateNum && hHmm.getParent() == parent) {
				return hHmm.getHmm();
			}
		}
		System.out.println("HMM not found!");
		return null;
	}

	public boolean hasDeeper(int depth, int parent) {
		if (depth + 1 > this.hmmSet.size() - 1)
			return false;
		List<HierarchicalHmm> nextLevel = this.hmmSet.get(depth+1).getHmmList();
		for (HierarchicalHmm hHmm : nextLevel) {
			if (hHmm.getParent() == parent)
				return true;
		}
		return false;
	}
}