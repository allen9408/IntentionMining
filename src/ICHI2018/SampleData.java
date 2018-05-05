package ICHI2018;

import HMM.Observation.ObservationDiscrete;

/*
 * This is an example from the paper "integrating Machine Learning and Workflow Management to .."
 */
public enum SampleData {
		start,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,end;
	public ObservationDiscrete<SampleData> observation() {
		return new ObservationDiscrete<SampleData>(this);
	}
}
