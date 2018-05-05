package ICHI2018;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import HMM.DataType.Hmm;
import HMM.Observation.ObservationDiscrete;
import HMM.PDF.OpdfDiscrete;

public class BuildInitialHMM {
	public ArrayList<SampleData> Observations = new ArrayList<SampleData>();
	public Hashtable<SampleData, Integer> activityIndex = new Hashtable<SampleData, Integer>();
	public Hashtable<String, ArrayList<Integer>> activityStringIndex = new Hashtable<String, ArrayList<Integer>>();
	public Hmm<ObservationDiscrete<SampleData>> hmm;

	public BuildInitialHMM(List<List<ObservationDiscrete<SampleData>>> sequences) {
		/* Construct an Observation List */
		for (List<ObservationDiscrete<SampleData>> seq : sequences) {
			for (ObservationDiscrete<SampleData> obj : seq) {
				if (!Observations.contains(obj.value)) {
					Observations.add(obj.value);
				}
			}
		}

		/* Construct HMM */
		hmm = new Hmm<ObservationDiscrete<SampleData>>(Observations.size());

		/* initialize activity index */
		int index = 0;
		for (SampleData obs : Observations) {
			activityIndex.put(obs, index);
			ArrayList<Integer> indexes = new ArrayList<Integer>();
			indexes.add(index);
			activityStringIndex.put(obs.toString(), indexes);
			index++;
		}

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
		for (int i = 0; i < Observations.size(); i++) {
			double[] stateObservationProb = new double[SampleData.values().length];
			for (int j = 0; j < SampleData.values().length; j++) {
				if (Observations.get(i) == SampleData.values()[j]) {
					stateObservationProb[j] = 1;
				}
			}
			hmm.setOpdf(i, new OpdfDiscrete<SampleData>(SampleData.class, stateObservationProb));
		}

		/* Set transition probability for states */
		hmm.setA(getTransformMatrix(sequences, Observations));
	}

	public double[][] getTransformMatrix(List<List<ObservationDiscrete<SampleData>>> sequences,
			ArrayList<SampleData> Observations) {
		int[][] matrix;

		// ~~~initialize the matrix
		if (0 != Observations.size()) {
			matrix = new int[Observations.size()][Observations.size()];
		} else {
			return null;
		}
		// ~~~Compute the Transform Matrix
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			SampleData previousActivity = SampleData.start;
			/* manually define the first activity */
			for (int i = 1; i < sequence.size(); i++) {
				SampleData currentActivity = sequence.get(i).value;
				int rowIndex = activityIndex.get(previousActivity);
				try {
					int columnIndex = activityIndex.get(currentActivity);
//					if(previousActivity.equals(SampleData.Airway_Verbalized)){
//						System.out.println(rowIndex+"--> " + columnIndex);
//					}
					matrix[rowIndex][columnIndex]++;
					previousActivity = currentActivity;
				} catch (NullPointerException e) {
					System.out.println(currentActivity);
				}
			}
		}
		return PrintandUpdateTransformMatrix(matrix, Observations);
	}

	private double[][] PrintandUpdateTransformMatrix(int[][] matrix, ArrayList<SampleData> Observations) {
		double[][] updatedMatrix = new double[Observations.size()][Observations.size()];
//		 System.out.println("The most general model");
		 for (SampleData obs : Observations) {
//			 System.out.print("\t" + obs + ",");
		 }

		 System.out.println();
		for (int i = 0; i < Observations.size(); i++) {
			for (int j = 0; j < Observations.size(); j++) {
				int sum = 0;
				for (int k = 0; k < Observations.size(); k++) {
					sum = sum + matrix[i][k];
				}
				if (0 != sum) {
					updatedMatrix[i][j] = (double) matrix[i][j] / sum;
//					 System.out.print("\t" + updatedMatrix[i][j]);
				} else {
//					 System.out.print("\t" + (double) matrix[i][j]);
				}
			}
//			 System.out.println();
		}
		return updatedMatrix;
	}

}