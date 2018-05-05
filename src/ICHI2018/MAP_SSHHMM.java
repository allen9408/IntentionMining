package ICHI2018;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/*
 * This method is based on the maximum likelihood
 * User can predetermine the state number or just let the program run to a point 
 * that furthur splitting won't increaset the probability any more
 * 
 * Written by Sen Yang 

 */

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;


import HMM.DataType.Hmm;
import HMM.Observation.Observation;
import HMM.Observation.ObservationDiscrete;
import HMM.PDF.OpdfDiscrete;
import HMM.PDF.OpdfDiscreteFactory;
import HMM.draw.GenericHmmDrawerDot;
import HMM.learn.BaumWelchScaledLearner;
import HMM.PathFindingAlgorithm.ViterbiCalculator;

public class MAP_SSHHMM<O extends Observation> {
	/************************ configuration parameters ***********************/

	int targetSize = Integer.MAX_VALUE;
	private static String fileName = "data/synthetic.csv";
	private static double printTransitionThreshold = 0;
	private Hmm<ObservationDiscrete<SampleData>> finalmodel;
	private static Hmm<ObservationDiscrete<SampleData>> iniHmm;
	private double lv;
	private static List<Integer> deepLevelActivityNumber = new ArrayList<Integer>();
	private static HierarchicalHmmSet hmmSet = new HierarchicalHmmSet();
	private static double lambda = 1;
	private static double lambda_1 = 1;


	private static Hashtable<Integer, AttributeDataSequences> attributeSet = new Hashtable<Integer, AttributeDataSequences>();
	private static List<StateAttributeSequence> stateAttributeSet;
//	private static boolean withAttribute = false;

	/************************ configuration parameters ***********************/
	public MAP_SSHHMM(String fileName) {
		this.fileName = fileName;
	}

	public MAP_SSHHMM() {
	}

	static public void main(String[] argv) throws java.io.IOException, CloneNotSupportedException {
		/* Input Data */
		attributeSet = getAttributeFromFile();
		SampleDataSequences sampleData = new SampleDataSequences(fileName);
		List<List<ObservationDiscrete<SampleData>>> sequences = sampleData.getSampleDataSequences();

		List<String> ID_set = new ArrayList<>(sampleData.get_IDlist());
		MAP_SSHHMM STACTsss = new MAP_SSHHMM();

		NoCrossValidation(sequences, ID_set, STACTsss);
		
		/* print stats */
		printingStats();
	}


	public static void NoCrossValidation(List<List<ObservationDiscrete<SampleData>>> sequences, List<String> ID_set, MAP_SSHHMM STACTsss) throws CloneNotSupportedException, IOException {
		lambda_1 = 1.0;
		lambda = 1.0;
		hmmSet.clear();
		STACTsss.STACTSSS(sequences, ID_set, attributeSet);
		// Output Average state numbers
		int totalSize = 0;
		for (Integer size : deepLevelActivityNumber) {
			totalSize += size;
		}
		int stateNumber = 0;
		int tranistionNumber = 0;

		for (int i = 0; i < hmmSet.size(); i++) {
			stateNumber += hmmSet.levelSize(i);
			for (Hmm<ObservationDiscrete<SampleData>> hmm : hmmSet.getLevelList(i)) {
				double[][] a = hmm.getA();
				for (int k = 0; k < hmm.nbStates(); k++) {
					for (int l = 0; l < hmm.nbStates(); l++) {
						if (a[k][l] != 0) {
							tranistionNumber++;
						}
					}
				}
			}
		}
		System.out.println(lambda + "," + hmmSet.size() + "," + stateNumber + "," + 1.0 * totalSize / deepLevelActivityNumber.size() + "," + tranistionNumber);
	}

	public static void printingStats() {
		// Output Average state numbers
		int totalSize = 0;
		for (Integer size : deepLevelActivityNumber) {
			totalSize += size;
		}
		int stateNumber = 0;
		int tranistionNumber = 0;
		System.out.println("Average number of activities in final level states: "
				+ 1.0 * totalSize / deepLevelActivityNumber.size());
		System.out.println("Depth = " + hmmSet.size());
		for (int i = 0; i < hmmSet.size(); i++) {
			stateNumber += hmmSet.levelSize(i);
			for (Hmm<ObservationDiscrete<SampleData>> hmm : hmmSet.getLevelList(i)) {
				double[][] a = hmm.getA();
				for (int k = 0; k < hmm.nbStates(); k++) {
					for (int l = 0; l < hmm.nbStates(); l++) {
						if (a[k][l] != 0) {
							tranistionNumber++;
						}
					}
				}
			}
		}
		System.out.println("Number of states: " + stateNumber);
		System.out.println("Average Number of states in each layer: " + stateNumber / hmmSet.size());
		System.out.println("Number of trainstions: " + tranistionNumber);
	}


	public static Hashtable<Integer, AttributeDataSequences> getAttributeFromFile() {
		//		String attributesFilename = "dataForICHI/attributes_05.15.2017.csv";
		String attributesFilename = "data/attributes.csv";
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(attributesFilename));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		try {
			reader = new BufferedReader(new FileReader(attributesFilename));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		Hashtable<Integer, AttributeDataSequences> attributeSet = new Hashtable<Integer, AttributeDataSequences>();
		try {
			String line = null;
			reader.readLine(); // the title
			while ((line = reader.readLine()) != null) {
				String[] record = line.split(",");
				AttributeDataSequences tmp = new AttributeDataSequences(Integer.parseInt(record[0]));
				for (int i = 0; i < 17; i++) {
					tmp.setAttribute(i, Integer.parseInt(record[i + 1]));
				}

				// Add attribute sequence to hash table
				attributeSet.put(Integer.parseInt(record[0]), tmp);

			}
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error in CsvFileReader !!!");
		}
		return attributeSet;
	}

	public <O extends Observation> void STACTSSS(List<List<ObservationDiscrete<SampleData>>> sequences,
												 List<String> ID_set, Hashtable<Integer, AttributeDataSequences> attributeSet)
			throws CloneNotSupportedException, IOException {

		/* Build Initial HMM */
		BuildInitialHMM_OneNodeAsInitial iniHMMbuilder = new BuildInitialHMM_OneNodeAsInitial(sequences);
		iniHmm = iniHMMbuilder.hmm.clone();

		/*
		 * SSS algorithm with Cross-entropy/Maximum Likelihood, with a specific
		 * tartget size
		 */
		finalmodel = SSSwithBIC(iniHmm, targetSize, sequences);

		/* Calculate NLL */
		printLogLikelihood(finalmodel, sequences);


		// Calculate Attribute probability in Transitions
		StateAttributeSequence[][] transitionAttributeSet = getAttributeProbs(finalmodel, sequences, ID_set,
				attributeSet);

		// Output transition attribute probabilities
		double[][] A = finalmodel.getA();


		hmmSet.addHmm(finalmodel, 0, 0, 0);
		// Output State Attribute Matrix
		String attributeOutfile = "SubStates/Attribute_0.md";
		stateAttributeSet = getStateAttributeProbs(finalmodel, sequences, ID_set, attributeSet);
		BufferedWriter attributeOut = new BufferedWriter(new FileWriter(attributeOutfile));
		for (int l = 0; l < finalmodel.nbStates(); l++) {
			attributeOut.write("##----State " + l + "----\n");
			StateAttributeSequence tmp = stateAttributeSet.get(l);
			for (int i = 0; i < 17; i++) {
				attributeOut.write("####Attribute " + i + "\n");
				attributeOut.write("| Value | Probability |\n");
				attributeOut.write("| ----- | ----------- |\n");

				ArrayList<Integer> attributes = tmp.getAttributeValueSequence(i);

				for (Integer x : attributes) {
					attributeOut.write("|" + x + " |" + tmp.getAttributeProbability(i, x) / tmp.size(i) + "|\n");
				}
			}
		}
		// Calculate Attribute probability in Transitions
		transitionAttributeSet = getAttributeProbs(finalmodel, sequences, ID_set, attributeSet);
		A = finalmodel.getA();

		for (int i = 0; i < finalmodel.nbStates(); i++) {
			for (int j = 0; j < finalmodel.nbStates(); j++) {
				if (A[i][j] != 0) {
					if (transitionAttributeSet[i][j].size(0) != 0) {
						attributeOut.write("##----Transition from " + i + " To " + j + "----\n");
						for (int k = 0; k < 17; k++) {
							attributeOut.write("####Attribute " + k + "\n");
							attributeOut.write("| Value | Probability |\n");
							attributeOut.write("| ----- | ----------- |\n");
							ArrayList<Integer> attributes = transitionAttributeSet[i][j].getAttributeValueSequence(k);
							for (Integer x : attributes) {
								attributeOut.write(
										"|" + x + " |" + transitionAttributeSet[i][j].getAttributeProbability(k, x)
												/ transitionAttributeSet[i][j].size(k) + "|\n");
							}
							attributeOut.write("");
						}
					}
				}
			}
		}

		// recurssively split hierarchical model
		hierarchicalSplitHMM(finalmodel, sequences, 1, 0, ID_set, "SubStates/SubState", attributeSet);
	}

	/**
	 * @param trainedHmm
	 * @param sequences
	 * @param ID_set
	 * @return
	 */
	public static List<StateAttributeSequence> getStateAttributeProbs(Hmm<ObservationDiscrete<SampleData>> trainedHmm,
																	  List<List<ObservationDiscrete<SampleData>>> sequences, List<String> ID_set,
																	  Hashtable<Integer, AttributeDataSequences> attributeSet) {
		// Path Finding -- Viterbi
		ArrayList<int[]> states = new ArrayList<int[]>();
		int index = 0;
		List<StateAttributeSequence> stateAttributeSet = new ArrayList<StateAttributeSequence>(trainedHmm.nbStates());
		for (int i = 0; i < trainedHmm.nbStates(); i++) {
			StateAttributeSequence tmp = new StateAttributeSequence(i);
			stateAttributeSet.add(tmp);
		}
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			ViterbiCalculator vc = new ViterbiCalculator(sequence, trainedHmm);
			states.add(vc.stateSequence());
		}
		for (int i = 0; i < states.size(); i++) {
			int[] state = states.get(i);
			int preState = state[0];
			int curState = preState;
			int idNum = Integer.parseInt(ID_set.get(i));
			// System.out.println("getID:" + idNum);
			AttributeDataSequences attributeSeq = attributeSet.get(idNum);
			if (attributeSeq == null) {
				System.out.println("NULL sequence");
			}
			StateAttributeSequence tmp = stateAttributeSet.get(curState);
			for (int k = 0; k < 17; k++) {
				// System.out.println(attributeSeq.getAttribute(k));
				tmp.addAttributeContent(k, attributeSeq.getAttribute(k));
			}
			for (int j = 1; j < state.length; j++) {
				// Label attributes to states
				// System.out.print(state[j] + ", ");
				curState = state[j];
				if (true) {
					tmp = stateAttributeSet.get(curState);
					for (int k = 0; k < 17; k++) {
						tmp.addAttributeContent(k, attributeSeq.getAttribute(k));
					}
				}
				preState = curState;
			}
		}
		return stateAttributeSet;
	}

	/**
	 * @param trainedHmm
	 * @param sequences
	 * @param ID_set
	 * @return
	 */
	public static StateAttributeSequence[][] getAttributeProbs(Hmm<ObservationDiscrete<SampleData>> trainedHmm,
															   List<List<ObservationDiscrete<SampleData>>> sequences, List<String> ID_set,
															   Hashtable<Integer, AttributeDataSequences> attributeSet) {
		// Path Finding -- Viterbi
		ArrayList<int[]> states = new ArrayList<int[]>();
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			ViterbiCalculator vc = new ViterbiCalculator(sequence, trainedHmm);
			states.add(vc.stateSequence());
		}
		// Calculate Attribute Probability Transition Matrix
		StateAttributeSequence[][] transitionAttributeSet = new StateAttributeSequence[trainedHmm.nbStates()][trainedHmm
				.nbStates()];
		for (int i = 0; i < trainedHmm.nbStates(); i++) {
			for (int j = 0; j < trainedHmm.nbStates(); j++) {
				transitionAttributeSet[i][j] = new StateAttributeSequence(j);
			}
		}
		for (int i = 0; i < states.size(); i++) {
			int[] state = states.get(i);
			int preState = state[0];
			int curState = preState;
			int idNum = Integer.parseInt(ID_set.get(i));
			// System.out.println("getID:" + idNum);
			AttributeDataSequences attributeSeq = attributeSet.get(idNum);
			if (attributeSeq == null) {
				System.out.println("NULL sequence");
			}

			for (int j = 1; j < state.length; j++) {
				// Label attributes to states
				curState = state[j];
				if (true) {
					for (int k = 0; k < 17; k++) {
						transitionAttributeSet[preState][curState].addAttributeContent(k, attributeSeq.getAttribute(k));
					}
				}
				preState = curState;
			}
		}
		return transitionAttributeSet;
	}

	/**
	 * @param trainedHmm
	 * @param sequences
	 * @return
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 */
	public static void hierarchicalSplitHMM(Hmm<ObservationDiscrete<SampleData>> trainedHmm,
											List<List<ObservationDiscrete<SampleData>>> sequences, int depth, int parent, List<String> ID_set,
											String subfile, Hashtable<Integer, AttributeDataSequences> attributeSet)
			throws IOException, CloneNotSupportedException {
		// Get state trace by viterbi
		ArrayList<int[]> states = new ArrayList<int[]>();
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			ViterbiCalculator vc = new ViterbiCalculator(sequence, trainedHmm);
			states.add(vc.stateSequence());
		}
		// No need to split start and end
		for (int i = 1; i < trainedHmm.nbStates(); i++) {
			if (i == 2) {
				// skip end
				continue;
			}
			String nextLevelFile = subfile + "_" + Integer.toString(i);
			String subStateFilename = nextLevelFile + ".csv";
			String subStateOutput = nextLevelFile + ".dot";
			String attributeOutfile = nextLevelFile + ".md";

			BufferedWriter substateOut = new BufferedWriter(new FileWriter(subStateFilename));
			// Create new sequence file for each state i
			String tabString = "";
			for (int j = 0; j < depth; j++) {
				tabString += "\t";
			}
			for (int j = 0; j < states.size(); j++) {
				int[] state = states.get(j);
				String id = ID_set.get(j);
				List<ObservationDiscrete<SampleData>> sequence = sequences.get(j);
				// System.out.print(id + ": ");
				for (int k = 0; k < state.length; k++) {
					if (state[k] == i) {
						substateOut.write(id + "," + sequence.get(k).value + "\n");
						// System.out.println(id + "," + sequence.get(k).value +
						// ",");
					}
				}
			}
			substateOut.close();
			// Using new acitivy traces to split sub-HMMs
			SampleDataSequences subsampleData = new SampleDataSequences(subStateFilename);
			List<List<ObservationDiscrete<SampleData>>> sub_sequences = subsampleData.getSampleDataSequences();
			List<String> sub_ID_set = new ArrayList<>(subsampleData.get_IDlist());
			MAP_SSHHMM STACTsss = new MAP_SSHHMM();
			// train sub model
			boolean goon = false;
			// System.out.println("ID Number:" + sub_ID_set.size());
			for (List<ObservationDiscrete<SampleData>> sequence : sub_sequences) {
				if (sequence.size() > 1) {
					goon = true;
				}
			}
			if (goon) {
				BuildInitialHMM_OneNodeAsInitial sub_iniHMMbuilder = new BuildInitialHMM_OneNodeAsInitial(
						sub_sequences);
				Hmm<ObservationDiscrete<SampleData>> sub_iniHmm = sub_iniHMMbuilder.hmm.clone();
				int targetSize = Integer.MAX_VALUE;
				Hmm<ObservationDiscrete<SampleData>> sub_finalmodel = SSSwithBIC(sub_iniHmm, targetSize, sub_sequences);
				// output training result
				List<ArrayList<Integer>> inState = new ArrayList<ArrayList<Integer>>(sub_finalmodel.nbStates());
				List<ArrayList<Double>> inStateProbs = new ArrayList<ArrayList<Double>>(sub_finalmodel.nbStates());
				for (int l = 0; l < sub_finalmodel.nbStates(); l++) {
					inState.add(new ArrayList<Integer>());
					inStateProbs.add(new ArrayList<Double>());
				}
				double[] divisionProb = new double[sub_finalmodel.nbStates()];
				for (int l = 0; l < sub_finalmodel.nbStates(); l++) {
					if (l == 2) {
						divisionProb[l] = 0;
						continue;
					}
					// System.out.println("-------Division: State: " + l + "
					// -------");
					// get division of each state
					double[] observationVector = new double[SampleData.values().length];
					int[] observationIndex = new int[SampleData.values().length];

					for (int j = 0; j < SampleData.values().length; j++) {
						observationVector[j] = sub_finalmodel.getOpdf(l)
								.probability(SampleData.values()[j].observation());
						observationIndex[j] = j;
					}

					// Sort observation vector
					for (int j = 0; j < observationVector.length - 1; j++) {
						for (int k = j + 1; k > 0; k--) {
							if (observationVector[k - 1] >= observationVector[k])
								break;
							double value = observationVector[k];
							int index = observationIndex[k];
							observationVector[k] = observationVector[k - 1];
							observationIndex[k] = observationIndex[k - 1];
							observationVector[k - 1] = value;
							observationIndex[k - 1] = index;
						}
					}
					List<Integer> division = new ArrayList<Integer>();
					int stateNum = 0;
					int divisionNum = 0;
					double divisionThreshold = 0.90;
					double probSum = 0.0;
					for (int j = 0; j < observationVector.length; j++) {
						if (observationVector[j] == 0.0)
							break;
						if (probSum < divisionThreshold) {
							inState.get(l).add(observationIndex[j]);
							inStateProbs.get(l).add(observationVector[j]);
							stateNum++;
							probSum += observationVector[j];
						} else {
							division.add(observationIndex[j]);
							divisionNum++;
						}
					}
				}

				hmmSet.addHmm(sub_finalmodel, depth, i, parent);
//				System.out.println("Add depth = " + depth + ", stateNum = " + i + ", parent = " + parent);
				(new GenericHmmDrawerDot()).writeSub(sub_finalmodel, subStateOutput, printTransitionThreshold, inState,
						inStateProbs);

				List<StateAttributeSequence> stateAttributeSet = getStateAttributeProbs(sub_finalmodel, sub_sequences,
						sub_ID_set, attributeSet);
				StateAttributeSequence[][] transitionAttributeSet = getAttributeProbs(sub_finalmodel, sub_sequences,
						sub_ID_set, attributeSet);
				BufferedWriter attributeOut = new BufferedWriter(new FileWriter(attributeOutfile));
				for (int l = 0; l < sub_finalmodel.nbStates(); l++) {
					attributeOut.write("##----State " + l + "----\n");
					StateAttributeSequence tmp = stateAttributeSet.get(l);
					for (int m = 0; m < 17; m++) {
						attributeOut.write("####Attribute " + m + "\n");
						attributeOut.write("| Value | Probability |\n");
						attributeOut.write("| ----- | ----------- |\n");

						ArrayList<Integer> attributes = tmp.getAttributeValueSequence(m);
						// System.out.println("Number of values: " +
						// tmp.getAttributeValueNumber(i));
						for (Integer x : attributes) {
							attributeOut
									.write("|" + x + " |" + tmp.getAttributeProbability(m, x) / tmp.size(m) + "|\n");
						}
					}
				}
				// Calculate Attribute probability in Transitions
				transitionAttributeSet = getAttributeProbs(sub_finalmodel, sequences, ID_set, attributeSet);
				double[][] A = sub_finalmodel.getA();

				for (int m = 0; m < sub_finalmodel.nbStates(); m++) {
					for (int n = 0; n < sub_finalmodel.nbStates(); n++) {
						if (A[m][n] != 0) {
							if (transitionAttributeSet[m][n].size(0) != 0) {
								attributeOut.write("##----Transition from " + m + " To " + n + "----\n");
								for (int o = 0; o < 17; o++) {
									attributeOut.write("####Attribute " + o + "\n");
									attributeOut.write("| Value | Probability |\n");
									attributeOut.write("| ----- | ----------- |\n");
									ArrayList<Integer> attributes = transitionAttributeSet[m][n]
											.getAttributeValueSequence(o);
									for (Integer x : attributes) {
										attributeOut.write("|" + x + " |"
												+ transitionAttributeSet[m][n].getAttributeProbability(o, x)
												/ transitionAttributeSet[m][n].size(o)
												+ "|\n");
									}
									attributeOut.write("");
								}
							}
						}
					}
				}
				attributeOut.close();

				if (sub_finalmodel.nbStates() > 3) {
					hierarchicalSplitHMM(sub_finalmodel, sub_sequences, depth + 1, i, sub_ID_set, nextLevelFile,
							attributeSet);
				} else {
					// System.out.println("#### Final " + nextLevelFile + "
					// activityNum: " + inState.get(1).size());
					deepLevelActivityNumber.add(inState.get(1).size());
				}
			}

		}
	}

	/**
	 * Calculate MAP
	 *
	 * @param trainedHmm
	 * @param sequences
	 * @return
	 */
	public static double calculateMapScore(Hmm<ObservationDiscrete<SampleData>> trainedHmm,
										   List<List<ObservationDiscrete<SampleData>>> sequences) {
		double score = 0.0;
		// calculate structure priors
		// q - number of states
		int q = trainedHmm.nbStates() + 1;
		// e - number of possible emmissions
		int e = SampleData.values().length + 1;
		// pStructure[i] - resulting structural contribution for state i
		double[] pStructure = new double[trainedHmm.nbStates()];
		// nT[i] - number of transitions from state i
		double[] nT = new double[trainedHmm.nbStates()];
		double[][] A = trainedHmm.getA();
		double nTSum = 0;
		for (int i = 0; i < A.length; i++) {
			nT[i] = 0;
			for (int j = 0; j < A.length; j++) {
				if (A[i][j] != 0 && !Double.isNaN(A[i][j])) {
					nT[i]++;
				}
			}
			nTSum += nT[i];
		}
		// nE[i] - numbers of emmission of state i
		double[] nE = new double[trainedHmm.nbStates()];
		double[] observationVector = new double[SampleData.values().length];
		double nESum = 0;
		for (int i = 0; i < trainedHmm.nbStates(); i++) {
			nE[i] = 0;
			if (i == 2)
				continue;
			nE[i] = 0;
			for (int j = 0; j < SampleData.values().length; j++) {
				observationVector[j] = trainedHmm.getOpdf(i).probability(SampleData.values()[j].observation());
				if (observationVector[j] != 0 && !Double.isNaN(observationVector[j])) {
					nE[i]++;
				}
			}
			nESum += nE[i];
		}
		double pt = nTSum / (q * trainedHmm.nbStates());
		double pe = nESum / (e * trainedHmm.nbStates());
		// Calculate structural priors

		for (int i = 1; i < trainedHmm.nbStates(); i++) {
			if (i == 2) {
				pStructure[i] = 0.0;
				continue;
			}

			pStructure[i] = Math.log(pt) * nT[i] + Math.log(1 - pt) * (q - nT[i])
					+ lambda_1 * (Math.log(pe) * nE[i] + Math.log(1 - pe) * (e - nE[i]));

		}

		// Calculate parameter priors
		double[] pParameter = new double[trainedHmm.nbStates()];
		// System.out.println("dirichletDist: " + dirichletDist);

		for (int i = 0; i < trainedHmm.nbStates(); i++) {
			if (i == 2) {
				pParameter[i] = 0;
				continue;
			}
			pParameter[i] = 0.0;
			double transitionProbs = 0.0;
			double transitionSum = 0.0;
			double emmisionProbs = 0.0;
			double emmisionSum = 0.0;
			for (int j = 0; j < trainedHmm.nbStates(); j++) {
				if (A[i][j] != 0 && !Double.isNaN(A[i][j])) {
					// System.out.println("A[i][j]: " + A[i][j]);
					transitionProbs += Math.log(A[i][j]);
					transitionSum += A[i][j];
				}
			}
			// System.out.println("State: " + i + "observationVector: ");
			for (int j = 0; j < SampleData.values().length; j++) {
				observationVector[j] = trainedHmm.getOpdf(i).probability(SampleData.values()[j].observation());
				if (observationVector[j] != 0 && !Double.isNaN(observationVector[j])) {
					// System.out.print(/*"observationVector[j]: " +
					// */observationVector[j] + ",");
					emmisionProbs += Math.log(observationVector[j]);
					emmisionSum += observationVector[j];
				}
			}
			// System.out.println("");
			pParameter[i] = transitionProbs + lambda_1 * emmisionProbs;

		}
		double dirichletDist = 0.0;
		for (int i = 0; i < trainedHmm.nbStates(); i++) {
			// System.out.println("pParameter: " + pParameter[i] + " pStructure:
			// " + pStructure[i]);
			// score += lambda_1 * pParameter[i] + pStructure[i];
			score += pParameter[i] + pStructure[i];
			if (nT[i] != 0 && nE[i] != 0) {
				// System.out.println("State: " + i + " " + "nT: " + nT[i] + "
				// nE: " + nE[i]);
				// System.out.print("State: " + i + " ");
				dirichletDist += getDirichletDistribution(2.0, (int) nT[i])
						+ lambda_1 * getDirichletDistribution(2.0, (int) nE[i]);
			}
		}
		score -= dirichletDist;


		// Calculate log(P(X|M))
		double loglikelihood = 0;
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			loglikelihood += trainedHmm.lnProbability(sequence);
		}

		// Calculate log(P(M)) + log(P(X|M))

		score = lambda * score + loglikelihood;
		return score;
	}

	/**
	 * Calculate Dirichlet Distribution
	 *
	 * @param alpha
	 * @param size
	 * @return
	 */
	public static double getDirichletDistribution(double alpha, int size) {
		double pro = 0.0;
		double numerator = 0.0;
		double denominator = 1.0;
		for (int i = 0; i < size; i++) {
			double tmp = (alpha - 0.5) * Math.log(alpha + 4.5) - (alpha + 4.5);
			double ser = 1.0 + 76.18009173 / (alpha + 0) - 86.50532033 / (alpha + 1) + 24.01409822 / (alpha + 2)
					- 1.231739516 / (alpha + 3) + 0.00120858003 / (alpha + 4) - 0.00000536382 / (alpha + 5);
			tmp += Math.log(ser * Math.sqrt(2 * Math.PI));
			numerator += tmp;
		}
		double tmp = (alpha * size - 0.5) * Math.log(alpha * size + 4.5) - (alpha * size + 4.5);
		double ser = 1.0 + 76.18009173 / (alpha * size + 0) - 86.50532033 / (alpha * size + 1)
				+ 24.01409822 / (alpha * size + 2) - 1.231739516 / (alpha * size + 3)
				+ 0.00120858003 / (alpha * size + 4) - 0.00000536382 / (alpha * size + 5);
		denominator = tmp + Math.log(ser * Math.sqrt(2 * Math.PI));
		pro = numerator - denominator;
		return pro;
	}

	/**
	 * Calculate MDL
	 *
	 * @param trainedHmm
	 * @param sequences
	 * @return
	 */
	public static double calculateMDL(Hmm<ObservationDiscrete<SampleData>> trainedHmm,
									  List<List<ObservationDiscrete<SampleData>>> sequences) {
		double activityNum = 0;
		ArrayList<ObservationDiscrete<SampleData>> traces = new ArrayList<ObservationDiscrete<SampleData>>();
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			int tempSum = 0;
			for (ObservationDiscrete<SampleData> activity : sequence) {
				tempSum++;
				if (!traces.contains(activity)) {
					traces.add(activity);
				}
			}
			activityNum += tempSum;
		}
		activityNum = Math.sqrt(activityNum);
		double MDLscoreTemp = 0;
		boolean ProIsValid = true;
		double Pro = 0;
		// as default, we set d =3
		int d = 3;
		int idx = 0;
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			Pro -= trainedHmm.lnProbability(sequence) / Math.log(2);
			if (Double.isNaN(Pro)) {
				ProIsValid = false;
				break;
			} else {
				Pro -= trainedHmm.lnProbability(sequence) / Math.log(2);
			}
		}

		if (ProIsValid) {
			double L = 0;
			int nA = traces.size();
			int nS = trainedHmm.nbStates();
			double nT = 0;
			double nB = 0;
			double[][] A = trainedHmm.getA();
			for (int i = 0; i < A.length; i++) {
				int tempSum = 0;
				for (int j = 0; j < A.length; j++) {
					if (A[i][j] > 0.01) {
						// nT++;
						tempSum++;
					}
				}
				nT += tempSum * tempSum;
			}
			nT = Math.sqrt(nT) / 4;
			for (int j = 0; j < trainedHmm.nbStates(); j++) {
				// System.out.print("State " + j + ": ");
				double[] observationVector = new double[SampleData.values().length];
				int tempSum = 0;
				for (int i = 0; i < SampleData.values().length; i++) {
					observationVector[i] = trainedHmm.getOpdf(j).probability(SampleData.values()[i].observation());
					// System.out.print(observationVector[i] + ", ");
					if (observationVector[i] != 0) {
						tempSum++;
					}
				}
				nB += tempSum * tempSum;
			}
			nB = Math.sqrt(nB);

			L = (activityNum + nB + nS) * Math.sqrt(nS) * (0 + Math.log(d + nA + 1)) / Math.log(2);
			System.out.println("nB: " + nB + "First part: " + (activityNum + nB) * Math.sqrt(nS) + " second part: "
					+ (0 + Math.log(d + nA + 1)) / Math.log(2));

			MDLscoreTemp = L + Pro; // not sure if it is necessary to multiply
			// the sequence number (I think not, since
			// the MDL is penalize the model)

			System.out.println("L: " + L + " pro: " + Pro);
		}
		return MDLscoreTemp;
	}

	private static double Lstar(int n) {
		double Lstar = Math.log(2.865) / Math.log(2);
		double N = n;
		while (Math.log(N) / Math.log(2) > 0) {
			Lstar += Math.log(N) / Math.log(2);
			N = Math.log(N) / Math.log(2);
		}
		return Lstar;
	}

	private static double binomial(int n, int k) {
		if (k > n - k)
			k = n - k;

		double logb = 0;
		for (int i = 1, m = n; i <= k; i++, m--)
			logb = logb + Math.log(m / i);
		return logb;
	}

	/**
	 * Calculate BIC Score
	 *
	 * @param trainedHmm
	 * @param sequences
	 * @return
	 */
	public static double calculateBICscore(Hmm<ObservationDiscrete<SampleData>> trainedHmm,
										   List<List<ObservationDiscrete<SampleData>>> sequences) {
		double loglikelihood = 0;
		int index = 0;
		int activityNum = 0;
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			loglikelihood += trainedHmm.lnProbability(sequence);
			// System.out.println("sequence " + index + " :" +
			// trainedHmm.lnProbability(sequence));
			index++;
			for (ObservationDiscrete<SampleData> activity : sequence) {
				activityNum++;
			}
		}

		/* Calculate Number of Free Parameters */
		double freeParameterNb = 0;
		double[] tempSum = new double[trainedHmm.nbStates()];
		/* free parameters in transition matrix A */
		double[][] A = trainedHmm.getA();
		for (int i = 0; i < A.length; i++) {
			// int tempSum = 0;
			for (int j = 0; j < A.length; j++) {
				if (A[i][j] != 0) {
					tempSum[i]++;
				}
			}
			// freeParameterNb += tempSum*tempSum;
		}
		/* free parameters in observation matrix B */
		for (int j = 0; j < trainedHmm.nbStates(); j++) {
			// System.out.print("State " + j + ": ");
			double[] observationVector = new double[SampleData.values().length];
			// int tempSum = 0;
			for (int i = 0; i < SampleData.values().length; i++) {
				observationVector[i] = trainedHmm.getOpdf(j).probability(SampleData.values()[i].observation());
				// System.out.print(observationVector[i] + ", ");
				if (observationVector[i] != 0) {
					tempSum[j]++;
				}
			}
			// freeParameterNb += tempSum*tempSum;
		}
		for (int i = 0; i < trainedHmm.nbStates(); i++) {
			freeParameterNb += tempSum[i] * tempSum[i] * trainedHmm.nbStates();
		}
		freeParameterNb = Math.sqrt(freeParameterNb);
		System.out.println("loglikelihood: " + loglikelihood + "  (freeParameterNb)*logT: "
				+ (freeParameterNb) * Math.log(activityNum));
		loglikelihood = loglikelihood - (freeParameterNb) * Math.log(activityNum);
		return loglikelihood;
	}

	public static void printLogLikelihood(Hmm<ObservationDiscrete<SampleData>> Hmm,
										  List<List<ObservationDiscrete<SampleData>>> sequences) {
		double sum = 0;
		int index = 0;
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			double prob = Hmm.lnProbability(sequence);
			// System.out.println("sequence " + index + " :" +
			// prob);
			sum += prob;
			index++;
		}
	}

	public Hmm<ObservationDiscrete<SampleData>> get_finalmodel() {
		return finalmodel;
	}

	public Hmm<ObservationDiscrete<SampleData>> get_inimodel() {
		return iniHmm;
	}

	public double getLv() {
		return lv;
	}

	/**
	 * SSS method
	 *
	 * @param learntHmm
	 * @return
	 */
	private static Hmm<ObservationDiscrete<SampleData>> SSSwithBIC(Hmm<ObservationDiscrete<SampleData>> learntHmm,
																   int targetSize, List<List<ObservationDiscrete<SampleData>>> sequences) {
		int iteration = 1;
		while (learntHmm.nbStates() < targetSize) {
			/*
			 * two parameters to store the best results in the current iteration
			 */
			Hmm<ObservationDiscrete<SampleData>> bestHmm = null;
			double bestLoglikelihood = Double.NEGATIVE_INFINITY;
			int bestSplitIndex = -1;

			/* the start state and end state should be excluded */
			for (int i = 1; i < learntHmm.nbStates(); i++) {
				/*
				 * state i and all its incoming and outgoing transitions are
				 * removed and new states i and learnHmm.nbstates() are inserted
				 * give random values to such transitions train the parameters
				 * with all the observations
				 */
				Hmm<ObservationDiscrete<SampleData>> copyHmm = new Hmm<ObservationDiscrete<SampleData>>(
						learntHmm.nbStates() + 1, new OpdfDiscreteFactory<SampleData>(SampleData.class));
				/* split state */
				iniSplitHmm(copyHmm, learntHmm, i);
				// (new GenericHmmDrawerDot()).print(copyHmm);

				/* train split model */
				long startTime = System.currentTimeMillis();
				Hmm<ObservationDiscrete<SampleData>> trainedHmm = trainSplitHmm(copyHmm, sequences, i);
				long trainTime = System.currentTimeMillis();
				// System.out.println("trainTime: " + (double)(trainTime -
				// startTime)/1000 + "s");

				/* calculate the log-likelihood */
				/* multiple observation sequence [Rabiner, p17] */
				// double loglikelihood = calculateBICscore(trainedHmm,
				// sequences);
				// double loglikelihood = calculateMDL(trainedHmm, sequences);
				double loglikelihood = calculateMapScore(trainedHmm, sequences);

				// System.out.println("state to split:" + i +" BIC: " +
				// loglikelihood);

				/*
				 * compare to the previous best results, update the best results
				 * if needed
				 */
				if (loglikelihood > bestLoglikelihood) {
					bestLoglikelihood = loglikelihood;
					bestSplitIndex = i;
					try {
						bestHmm = trainedHmm.clone();
					} catch (CloneNotSupportedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			// System.out.println("BestSplitIndex: " + bestSplitIndex + " BIC: "
			// + bestLoglikelihood);

			/* calculate the log-likelihood with the old model */

			double oldLoglikelihood = calculateMapScore(learntHmm, sequences);
			if (bestLoglikelihood / oldLoglikelihood > 1) {
				/*
				 * break the iteration when the log-likelihood of new model is
				 * smaller than its last-step model
				 */
				break;
			} else {
				/* update the current model with the new model */
				try {
					learntHmm = bestHmm.clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			iteration++;
			// System.out.println("OldBIC: " + oldLoglikelihood + " BestBIC: " +
			// bestLoglikelihood);
		}

		// smoothing
		learntHmm = smoothing(learntHmm, sequences);

		// division
		double[] divisionProb = getDivision(learntHmm);

		return learntHmm;
	}

	// Laplace smoothing
	public static Hmm<ObservationDiscrete<SampleData>> smoothing(Hmm<ObservationDiscrete<SampleData>> trainedHmm,
																 List<List<ObservationDiscrete<SampleData>>> sequences) {
		for (int i = 1; i < trainedHmm.nbStates(); i++) {
			if (i == 2) {
				continue;
			}
			double[] observationVector = new double[SampleData.values().length];
			int n = sequences.size();
			// int n = SampleData.values().length;
			int ts = trainedHmm.nbStates();
			double v = 0;
			for (int j = 0; j < SampleData.values().length; j++) {
				observationVector[j] = trainedHmm.getOpdf(i).probability(SampleData.values()[j].observation());
				if (observationVector[j] > 0) {
					v++;
				}
			}
			double p = 0.005 / (SampleData.values().length - 2 - v);
			double q = 0.005 / v;
			for (int j = 1; j < SampleData.values().length - 1; j++) {
				if (observationVector[j] > 0) {
					observationVector[j] = observationVector[j] * (1 - 0.005);
				} else {
					observationVector[j] += p;
				}
				// observationVector[j] = (n * observationVector[j] + 1) / (n +
				// SampleData.values().length - 2);

				// System.out.println("Observation Vector: " +
				// observationVector[j]);
			}
			trainedHmm.setOpdf(i, new OpdfDiscrete<SampleData>(SampleData.class, observationVector));
		}
		return trainedHmm;
	}

	public static double[] getDivision(Hmm<ObservationDiscrete<SampleData>> trainedHmm) {
		double[] divisionProb = new double[trainedHmm.nbStates()];
		for (int i = 0; i < trainedHmm.nbStates(); i++) {
			if (i == 2) {
				divisionProb[i] = 0;
				continue;
			}
			// System.out.println("-------Division: State: " + i + " -------");
			// get division of each state
			double[] observationVector = new double[SampleData.values().length];
			int[] observationIndex = new int[SampleData.values().length];

			for (int j = 0; j < SampleData.values().length; j++) {
				observationVector[j] = trainedHmm.getOpdf(i).probability(SampleData.values()[j].observation());
				observationIndex[j] = j;
			}

			// Sort observation vector
			for (int j = 0; j < observationVector.length - 1; j++) {
				for (int k = j + 1; k > 0; k--) {
					if (observationVector[k - 1] >= observationVector[k])
						break;
					double value = observationVector[k];
					int index = observationIndex[k];
					observationVector[k] = observationVector[k - 1];
					observationIndex[k] = observationIndex[k - 1];
					observationVector[k - 1] = value;
					observationIndex[k - 1] = index;
				}
			}
			List<Integer> state = new ArrayList<Integer>();
			List<Integer> division = new ArrayList<Integer>();
			int stateNum = 0;
			int divisionNum = 0;
			double divisionThreshold = 0.90;
			double probSum = 0.0;
			for (int j = 0; j < observationVector.length; j++) {
				if (observationVector[j] == 0.0)
					break;
				if (probSum < divisionThreshold) {
					state.add(observationIndex[j]);
					stateNum++;
					probSum += observationVector[j];
				} else {
					division.add(observationIndex[j]);
					divisionNum++;
				}
			}

			divisionProb[i] = 1 - probSum;
		}
		return divisionProb;
	}

	public static double getBetaDistribution(double alpha, double beta) {
		double a = alpha + beta;
		double b;
		if (Math.min(alpha, beta) <= 1) {
			b = Math.max(1 / alpha, 1 / beta);
		} else {
			b = Math.sqrt((a - 2) / (2 * alpha * beta - a));
		}
		double c = alpha + 1 / b;

		double W = 0;
		boolean reject = true;
		for (reject = true; reject; ) {
			double U1 = Math.random();
			double U2 = Math.random();
			double V = b * Math.log(U1 / (1 - U1));
			W = alpha * Math.exp(V);
			reject = (a * Math.log(a / (beta + W)) + c * V - Math.log(4)) < Math.log(U1 * U1 * U2);
		}
		return (W / (beta + W));
	}

	/**
	 * ]
	 *
	 * @param copyHmm
	 * @param learntHmm
	 * @param splitStateIndex
	 */
	public static void iniSplitHmm(Hmm<ObservationDiscrete<SampleData>> copyHmm,
								   Hmm<ObservationDiscrete<SampleData>> learntHmm, int splitStateIndex) {
		/* set Pi */
		copyHmm.setPi(0, 1);
		/* set the observation matrix B */
		for (int j = 0; j < learntHmm.nbStates(); j++) {
			/*
			 * set the observation of the newly added state equal to the ith
			 * state in learntHMM
			 */
			if (j == splitStateIndex) {
				/* the splitting-state observation distribution */
				double[] observationVector = new double[SampleData.values().length];
				double[] splittingStateObservationVector = new double[SampleData.values().length];
				double[] newSplitStateObservationVector = new double[SampleData.values().length];
				for (int i = 1; i < SampleData.values().length
						- 1; i++) { /*
									 * i!=0 && i!=SampleData.values().length-1
									 */
					observationVector[i] = learntHmm.getOpdf(splitStateIndex)
							.probability(SampleData.values()[i].observation());
					// splittingStateObservationVector[i] =
					// Math.abs(Math.random() - 0.5) * observationVector[i];
					splittingStateObservationVector[i] = getBetaDistribution(2.0, 2.0) * observationVector[i];
					newSplitStateObservationVector[i] = observationVector[i] - splittingStateObservationVector[i];
				}
				/* set the observation distribution of the splitting state */
				copyHmm.setOpdf(splitStateIndex,
						new OpdfDiscrete<SampleData>(SampleData.class, splittingStateObservationVector));
				/* set the distribution of the new split state */
				copyHmm.setOpdf(learntHmm.nbStates(),
						new OpdfDiscrete<SampleData>(SampleData.class, newSplitStateObservationVector));
			} else {
				double[] observationVector = new double[SampleData.values().length];
				for (int i = 0; i < SampleData.values().length; i++) {
					observationVector[i] = learntHmm.getOpdf(j).probability(SampleData.values()[i].observation());
				}
				copyHmm.setOpdf(j, new OpdfDiscrete<SampleData>(SampleData.class, observationVector));
			}
		}
		/* set the transition matrix */
		/* give random values to the transitions of the split state */
		double[][] transitionMatrix = new double[learntHmm.nbStates() + 1][learntHmm.nbStates() + 1];
		double[][] lastStepTransitionMatrix = new double[learntHmm.nbStates()][learntHmm.nbStates()];
		lastStepTransitionMatrix = learntHmm.getA();
		/* fill and perturb the transition matrix */
		fillAndPerturbTransitionMatrix(transitionMatrix, lastStepTransitionMatrix, splitStateIndex);
		// System.out.println("Split transition matrix:");
		// printTransitionMatrix(transitionMatrix);
		copyHmm.setA(transitionMatrix);
	}

	/**
	 * @param transitionMatrix
	 * @param lastStepTransitionMatrix
	 * @param splitStateIndex
	 */
	public static void fillAndPerturbTransitionMatrix(double[][] transitionMatrix, double[][] lastStepTransitionMatrix,
													  int splitStateIndex) {

		double rowISum = 0;
		double rowLastSum = 0;
		for (int i = 0; i < lastStepTransitionMatrix.length; i++) {
			for (int j = 0; j < lastStepTransitionMatrix.length; j++) {
				/* perturb the split state */
				if (lastStepTransitionMatrix[i][j] != 0) {
					/* fill with a random number */
					/*
					 * ??????? do i need to keep the sum of transitions to 1 ???
					 */
					/* also give a random number to the new state */
					if (i == splitStateIndex) { // Out (two states (A and A')
						// --> multiple)
						double new_tran1 = getBetaDistribution(2.0, 2.0) * lastStepTransitionMatrix[i][j];
						// System.out.println("Beta: " + new_tran1);
						double new_tran2 = lastStepTransitionMatrix[i][j] - new_tran1;
						transitionMatrix[lastStepTransitionMatrix.length][j] = new_tran1;
						transitionMatrix[i][j] = new_tran2;
						rowLastSum += new_tran1;
						rowISum += new_tran2;
					} else if (j == splitStateIndex) {
						double new_tran = getBetaDistribution(2.0, 2.0) * lastStepTransitionMatrix[i][j];
						transitionMatrix[i][lastStepTransitionMatrix.length] = new_tran;
						transitionMatrix[i][j] = lastStepTransitionMatrix[i][j] - new_tran;
					} else {
						transitionMatrix[i][j] = lastStepTransitionMatrix[i][j];
					}
				}
			}
		}
		/* add self transition */
		transitionMatrix[splitStateIndex][lastStepTransitionMatrix.length] = lastStepTransitionMatrix[splitStateIndex][splitStateIndex]
				* getBetaDistribution(2.0, 2.0);
		transitionMatrix[lastStepTransitionMatrix.length][lastStepTransitionMatrix.length] = lastStepTransitionMatrix[splitStateIndex][splitStateIndex]
				- transitionMatrix[splitStateIndex][splitStateIndex];
		rowLastSum += transitionMatrix[lastStepTransitionMatrix.length][lastStepTransitionMatrix.length];
		/* normalized last row */
		// Normalize the splitted row and last row
		if (rowISum != 0 && rowLastSum != 0) {
			for (int j = 0; j <= lastStepTransitionMatrix.length; j++) {
				transitionMatrix[lastStepTransitionMatrix.length][j] = transitionMatrix[lastStepTransitionMatrix.length][j]
						/ rowLastSum;
				transitionMatrix[splitStateIndex][j] = transitionMatrix[splitStateIndex][j] / rowISum;
			}
		}
	}

	/**
	 * @param copyHmm
	 * @param sequences
	 * @return
	 */
	public static Hmm<ObservationDiscrete<SampleData>> trainSplitHmm(Hmm<ObservationDiscrete<SampleData>> copyHmm,
																	 List<List<ObservationDiscrete<SampleData>>> sequences, int curState) {
		BaumWelchScaledLearner bwl = new BaumWelchScaledLearner();
		Hmm<ObservationDiscrete<SampleData>> trainedHmm = bwl.learn(copyHmm, sequences);


		double[] s1 = new double[SampleData.values().length];
		double[] s2 = new double[SampleData.values().length];
		for (int i = 0; i < SampleData.values().length; i++) {
			s1[i] = trainedHmm.getOpdf(curState).probability(SampleData.values()[i].observation());
			s2[i] = trainedHmm.getOpdf(trainedHmm.nbStates() - 1).probability(SampleData.values()[i].observation());

		}

		double sum_1 = 0.0;
		double sum_2 = 0.0;

		// balance vectors
		double threshold = 0.00001;
		/*
		 * start and end should not be considered this may be the reason that
		 * triggers NaN
		 */
		for (int i = 1; i < SampleData.values().length; i++) {
			if (i == 2)
				continue;
			// System.out.println("check observation: " + i);
			if ((s1[i] < threshold) || (s2[i] < threshold)) {
				// System.out.println("rebalance observation: " + i);
				if (s1[i] < s2[i]) {
					s1[i] = 0.0;
				} else {
					s2[i] = 0.0;
				}
			}
			sum_1 += s1[i];
			sum_2 += s2[i];
		}

		// System.out.println("sum1: " + sum_1 + "sum_2: " + sum_2);

		for (int i = 0; i < SampleData.values().length; i++) {
			s1[i] /= sum_1;
			s2[i] /= sum_2;
			// System.out.println("Observation - rebalance " + i + ": " + s1[i]
			// + "\t"+s2[i]);
		}

		// return balanced vactors to trained HMM
		trainedHmm.setOpdf(curState, new OpdfDiscrete<SampleData>(SampleData.class, s1));
		trainedHmm.setOpdf(trainedHmm.nbStates() - 1, new OpdfDiscrete<SampleData>(SampleData.class, s2));
		// (new
		// GenericHmmDrawerDot()).print(trainedHmm,printTransitionThreshold);

		Hmm<ObservationDiscrete<SampleData>> trainedHmm2 = bwl.learn(trainedHmm, sequences);
		// (new
		// GenericHmmDrawerDot()).print(trainedHmm2,printTransitionThreshold);

		return trainedHmm2;
	}
}