package HMM.IO;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import HMM.DataType.Hmm;
import HMM.Observation.Observation;
import HMM.Observation.ObservationDiscrete;
import HMM.PDF.OpdfDiscrete;
import ICHI2018.SampleData;

public class Import {
	public int nbstate;
	public double[][] transitionMatrix;
	public double[] pi;
	public ArrayList<Opdfs> ArrayOpdfs = new ArrayList<Opdfs>();
	public Import(String filename) throws FileNotFoundException {
		// TODO Auto-generated constructor stub
	      String  thisLine = null;
	          // open input stream test.txt for reading purpose.
	    	  FileInputStream is = new FileInputStream(filename);
	    	  InputStreamReader isr = new InputStreamReader(is);
	          BufferedReader br = new BufferedReader(isr);
	          try {
				while ((thisLine = br.readLine()) != null) {
					  if (thisLine.equals("nbstate")) {
						  if ((thisLine = br.readLine()) != null) {
				    		  nbstate = Integer.valueOf(thisLine);
				    	      transitionMatrix = new double[nbstate][nbstate];
				    	      pi = new double[nbstate];
						  }
					  }
					  else if (thisLine.equals("transitionMatrix")) {
						  int i = 0;
						  while ((thisLine = br.readLine()) != null && i < nbstate) {
							  String[] line = thisLine.split(" ");
							  int j = 0;
							  for (String s : line) {
								  if (!s.equals("")) {
									  transitionMatrix[i][j] = Double.valueOf(s);
									  j++;
								  }
							  }
							  i++;
						  }
					  }
					  else if (thisLine.equals("pi")) {
						  int i = 0;
						  while ((thisLine = br.readLine()) != null && i < nbstate) {
							  pi[i] = Double.valueOf(thisLine);
							  i++;
						  }
					  }
					  else if (thisLine.equals("opdfs")) {
//	        		  String first = "", second  = "";
//	        		  while ((first = br.readLine()) != null && 
//	        				  (second  = br.readLine())!= null) {
//	        			  String[] firstsplit = first.split(",");
//	        			  String[] secondsplit = second.split(",");
//	        			  Opdfs op = new Opdfs(firstsplit, secondsplit);
//	        			  opdfs.add(op);
//	        		  }
						  String first = "";
						  while ((first = br.readLine()) != null) {
							  String[] firstsplit = first.split(",");
							  double[] prob = new double[firstsplit.length];
							  int index = 0;
							  for(String split:firstsplit){
								  prob[index] = Double.valueOf(firstsplit[index]);
								  index++;
							  }
							  Opdfs op = new Opdfs(prob);
							  ArrayOpdfs.add(op);
						  }
					  }
				  }
			} catch (NumberFormatException | IOException |NullPointerException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	}
	
	public Hmm<ObservationDiscrete<SampleData>> getHMM(){
		Hmm<ObservationDiscrete<SampleData>> hmm= new Hmm<ObservationDiscrete<SampleData>>(nbstate);

		/*
		 * Initialize all the parameters for the HMM Pi, A, B
		 */

		/* Set Pi */
		/* use the first state, i.e. initial state, as the start state */
		hmm.setPi(0, 1);

		/* Set the observation matrix B */
		for (int i = 0; i < nbstate; i++) {
			hmm.setOpdf(i, new OpdfDiscrete<SampleData>(SampleData.class, ArrayOpdfs.get(i).getProb()));
		}
		System.out.println(SampleData.values().length);

		/* Set transition probability for states */
		hmm.setA(transitionMatrix);
		
		return hmm;
	}
	
	
	
	
}
