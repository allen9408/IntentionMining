package ICHI2018;


/*
 * This method is based on the maximum likelihood
 * User can predetermine the state number or just let the program run to a point 
 * that furthur splitting won't increaset the probability any more
 * 
 * Written by Sen Yang 

 */



import HMM.DataType.Hmm;
import HMM.Observation.ObservationDiscrete;
import ICHI2018.SampleData;


public class HierarchicalHmm {
	private int depth;
	private int stateNum;
	private int parent;
	private Hmm<ObservationDiscrete<SampleData>> Hmm;

	public HierarchicalHmm(Hmm<ObservationDiscrete<SampleData>> Hmm, int depth, int stateNum, int parent) {
		this.Hmm = Hmm;
		this.depth = depth;
		this.stateNum = stateNum;
		this.parent = parent;
	}

	public Hmm<ObservationDiscrete<SampleData>> getHmm() {
		return this.Hmm;
	}
	public int getDepth() {
		return this.depth;
	}
	public int getStateNum() {
		return this.stateNum;
	}
	public int getParent() {
		return this.parent;
	}
}