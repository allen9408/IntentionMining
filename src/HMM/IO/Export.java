package HMM.IO;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import HMM.DataType.Hmm;
import HMM.Observation.Observation;
import HMM.Observation.ObservationDiscrete;
import HMM.PDF.Opdf;

public class Export<O extends Observation> {

	public Export(Hmm<O> hmm, String modelType) {
		// TODO Auto-generated constructor stub
		/*
		 * four value that need to be exported to file
		 */
		try {
			// get current time
			DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss");
			Date date = new Date();
			String curr = dateFormat.format(date);
			// create a new file
			File file = new File(modelType + curr + ".txt");
			PrintWriter w = new PrintWriter(file, "UTF-8");
			// write nbstate
			w.println("nbstate");
			int nbstate = hmm.nbStates();
			w.println(nbstate);
			// write transition matrix
			w.println("transitionMatrix");
			double[][] transitionMatrix = hmm.getA();
			for (int i = 0; i < transitionMatrix.length; i++) {
				for (int j = 0; j < transitionMatrix[0].length; j++) {
					w.print(transitionMatrix[i][j]);
//					System.out.println(transitionMatrix[i][j]);
					w.print(" ");
				}
				w.println();
			}
			// write pi list
			w.println("pi");
			double[] pi = hmm.getPi();
			for (int i = 0; i < pi.length; i++) {
				w.println(pi[i]);
			}
			// write opdfs
			w.println("opdfs");
			ArrayList<Opdf<O>> opdfs = hmm.getOpdf();
			for (int i = 0; i < opdfs.size(); i++) {
				w.println(opdfs.get(i).toString());
			}
			w.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
}
