package ICHI2018;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import HMM.DataType.Hmm;
import HMM.Observation.ObservationDiscrete;
import HMM.PDF.OpdfDiscrete;
import HMM.draw.GenericHmmDrawerDot;

public class BuildInitialHMM_OneNodeAsInitial {
	public Hmm<ObservationDiscrete<SampleData>> hmm;
	public Hashtable<SampleData, Integer> ObservationAndFrequency = new Hashtable<SampleData, Integer>();
	private int activityNum = 0;

	public BuildInitialHMM_OneNodeAsInitial(List<List<ObservationDiscrete<SampleData>>> sequences) {
		for (List<ObservationDiscrete<SampleData>> seq : sequences) {
			for (ObservationDiscrete<SampleData> obj : seq) {
				activityNum++;
				if (!ObservationAndFrequency.keySet().contains(obj.value)) {
					ObservationAndFrequency.put(obj.value, 1);
				}else{
					ObservationAndFrequency.put(obj.value, ObservationAndFrequency.get(obj.value)+1);
				}
			}
			/* activity num without conunting start and end node */
			activityNum-=2;
		}
		
		/* Construct HMM */
		hmm = new Hmm<ObservationDiscrete<SampleData>>(3);

		buildMMInitHmm(sequences);
	}

	public void buildMMInitHmm(List<List<ObservationDiscrete<SampleData>>> sequences) {

		/*
		 * Initialize all the parameters for the HMM Pi, A, B
		 */

		/* Set Pi */
		/* use the first state, i.e. initial state, as the start state */
		hmm.setPi(0, 1);

		/* Set the observation matrix B */
		/* set the start node */
		double[] stateObservationProb_start = new double[SampleData.values().length];
		stateObservationProb_start[0] = 1;
		hmm.setOpdf(0, new OpdfDiscrete<SampleData>(SampleData.class, stateObservationProb_start));
		// System.out.println("SampleData.values(): " + SampleData.values().length);
		double[] stateObservationProb_firstNode = new double[SampleData.values().length];
		for(int i=1; i<SampleData.values().length-1; i++){
//			 System.out.println(i + "Iteration: Get: "+ SampleData.values()[i]);
			 if (ObservationAndFrequency.get(SampleData.values()[i]) != null) {
				stateObservationProb_firstNode[i] = 1.0*ObservationAndFrequency.get(SampleData.values()[i])/activityNum;
			 } else {
				 stateObservationProb_firstNode[i] = 0.0;
			 }
		}
		hmm.setOpdf(1, new OpdfDiscrete<SampleData>(SampleData.class, stateObservationProb_firstNode));	
		
		/* the end node */
		double[] stateObservationProb_end = new double[SampleData.values().length];
		stateObservationProb_end[SampleData.values().length-1] = 1;
		hmm.setOpdf(2, new OpdfDiscrete<SampleData>(SampleData.class, stateObservationProb_end));

		/* Set transition probability for states */
		hmm.setA(getTransitionMatrix(sequences));
		
		/* print initial matrix */
		// (new GenericHmmDrawerDot()).print(hmm,
		// 		 0);
	}

	public double[][] getTransitionMatrix(List<List<ObservationDiscrete<SampleData>>> sequences) {
		double[][]	matrix = new double[3][3];

		/* start to firstNode */
		matrix[0][1] = 1;
		
		/* firstNode to firstNode/End */
		int numEnd = 0;
		int numAct = 0;
		for(List<ObservationDiscrete<SampleData>> sequence:sequences){
			numEnd++;
			for(ObservationDiscrete<SampleData> activity:sequence){
				numAct++;
			}
		}
		/* firstNode selftransition */
		matrix[1][1] = (numAct-3*numEnd)*1.0/(numAct-2*numEnd);
		/* firstNode to End */
		matrix[1][2] = numEnd*1.0/(numAct-2*numEnd);
		


		return matrix;
	}

}
