package ICHI2018;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import HMM.Observation.ObservationDiscrete;

public class SampleDataSequences {
	private List<List<ObservationDiscrete<SampleData>>> sequences = new ArrayList<List<ObservationDiscrete<SampleData>>>();
	private ArrayList<ArrayList<ObservationDiscrete<SampleData>>> log = new ArrayList<ArrayList<ObservationDiscrete<SampleData>>>();
	int symbolTypeNumber;
	private Hashtable<SampleData, Integer> activityIndex = new Hashtable<SampleData, Integer>();
	private ArrayList<String> ID_set = new ArrayList<String>();
	private List<String> Activity_Set_Mapping = new ArrayList<String>(
		Arrays.asList("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o"));
	private String file_name = "";

	public SampleDataSequences(String fileName) {
		this.file_name = fileName;
		readCsvFile();
	}

	public SampleDataSequences() {
		file_name = "IntubationData_06.03.csv";
		readCsvFile();
	}

	public boolean readCsvFile() {

		boolean isSuccess = true;

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(file_name));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}

		ArrayList<ObservationDiscrete<SampleData>> trace = new ArrayList<ObservationDiscrete<SampleData>>();

		try {
			String line = null;
			reader.readLine(); // The title
			while ((line = reader.readLine()) != null) {
				String[] record = line.split(",");
				// record[0] = record.toLoerCase();
				if (ID_set.isEmpty()) { // Initialization
					trace.add(SampleData.start.observation());
					ID_set.add(record[0]);
//					System.out.print("ID :" + record[0] + "\t");
				} else {
					if (!ID_set.contains(record[0])) {
						if (!trace.isEmpty()) {
							// System.out.println(ID + ":");
							trace.add(SampleData.end.observation());
							ArrayList<ObservationDiscrete<SampleData>> temp_trace = (ArrayList<ObservationDiscrete<SampleData>>) trace
									.clone();
							log.add(temp_trace);
							sequences.add(temp_trace);
							// print_trace(trace);
							trace.clear();
							trace.add(SampleData.start.observation());
						}
						ID_set.add(record[0]);
//						System.out.print("ID :" + record[0] + "\t");
					}
				}
				record[1] = record[1].replace("(", "");
				record[1] = record[1].replace(")", "");
				if (!Activity_Set_Mapping.contains(record[1].toLowerCase())) {
					System.out.println("!!!ID: " + record[0] + "activity: " + record[1]);
				}
				// System.out.println(record[1]);
				int index = Activity_Set_Mapping.indexOf(record[1].toLowerCase());
				// System.out.println("Index" + index);
				SampleData event = SampleData.values()[index + 1];
				trace.add(event.observation());
			}
			trace.add(SampleData.end.observation());
			// print the last trace
			// print_trace(trace);
			log.add(trace);
			sequences.add(trace);
			// Last line
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error in CsvFileReader !!!");
			isSuccess = false;
			return false;
		} finally {
			if (isSuccess) {
				// print_IDset(ID_set);
				iniParameters();
				try {
					reader.close();
				} catch (IOException e) {
					System.out.println("Error while closing fileReader/csvFileParser !!!");
				}
			}
		}
		if (isSuccess) {

			// print_log();
			return true;
		}
		return false;
	}

	public void print_IDset(ArrayList<String> ID_set) {
		for (String ID : ID_set) {
			System.out.println(ID);
		}
	}

	public ArrayList<String> get_IDlist() {
		return this.ID_set;
	}

	public void print_log() {
		for (int i = 0; i < log.size(); i++) {
			print_trace(log.get(i));
		}
	}

	public void print_trace(List<ObservationDiscrete<SampleData>> trace) {
		for (int i = 0; i < trace.size(); i++) {
			System.out.print(trace.get(i)+", ");
		}
		System.out.println();
	}

	public ArrayList<ArrayList<ObservationDiscrete<SampleData>>> get_log() {
		return log;
	}

	public void iniParameters() {
		symbolTypeNumber = SampleData.values().length;

		/* initialize activity index */
		for (int i = 0; i < SampleData.values().length; i++) {
			activityIndex.put(SampleData.values()[i], i);
		}
//		System.out.println();
	}

	public List<List<ObservationDiscrete<SampleData>>> getSampleDataSequences() {
		return sequences;
	}
	public List<List<ObservationDiscrete<SampleData>>> getRecommendSequences() {
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			sequence.remove(sequence.size() - 1);
		}
		return sequences;
	}

	public int getSampleDataSymbolSize() {
		int size = 0;
		for (List<ObservationDiscrete<SampleData>> sequence : log) {
			/* only count the symbols, do not count the start and end state */
			size += sequence.size() - 2;
		}
		return size;
	}

	public double[][] getTransformMatrix() {
		int[][] matrix;
		// ~~~initialize the matrix
		if (0 != symbolTypeNumber) {
			matrix = new int[symbolTypeNumber][symbolTypeNumber];
		} else {
			return null;
		}
		// ~~~Compute the Transform Matrix
		for (List<ObservationDiscrete<SampleData>> sequence : sequences) {
			// print_trace(sequence);
			SampleData previousActivity = SampleData.start;
			SampleData currentActivity;

			/* manually define the first activity */
			for (int i = 1; i < sequence.size(); i++) {
				currentActivity = sequence.get(i).value;
				int rowIndex = activityIndex.get(previousActivity);
				int columnIndex = activityIndex.get(currentActivity);
				matrix[rowIndex][columnIndex] += 1;
				previousActivity = currentActivity;
			}
		}
		return PrintandUpdateTransformMatrix(matrix);
	}

	private double[][] PrintandUpdateTransformMatrix(int[][] matrix) {
		double[][] updatedMatrix = new double[symbolTypeNumber][symbolTypeNumber];
		for (int i = 0; i < SampleData.values().length; i++) {
			System.out.print("\t" + SampleData.values()[i]);
		}

		System.out.println();
		for (int i = 0; i < symbolTypeNumber; i++) {
			int sum = 0;
			for (int k = 0; k < symbolTypeNumber; k++) {
				sum = sum + matrix[i][k];
			}
			for (int j = 0; j < symbolTypeNumber; j++) {

				if (0 != sum) {
					updatedMatrix[i][j] = (double) matrix[i][j] / sum;
					System.out.print("\t" + updatedMatrix[i][j]);
				} else {
					System.out.print("\t" + (double) matrix[i][j]);
				}
			}
			System.out.println();
		}
		return updatedMatrix;
	}

}
