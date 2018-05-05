/* jahmm package - v0.6.1 */

/*
  *  Copyright (c) 2004-2006, Jean-Marc Francois.
 *
 *  This file is part of Jahmm.
 *  Jahmm is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Jahmm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Jahmm; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

package HMM.draw;

import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import HMM.DataType.Hmm;
import HMM.Observation.Observation;
import HMM.PDF.Opdf;
import ICHI2018.SampleData;

/**
 * An HMM to <i>dot</i> file converter. See
 * <url>http://www.research.att.com/sw/tools/graphviz/</url> for more
 * information on the <i>dot</i> tool.
 * <p>
 * The command <tt>dot -Tps -o &lt;outputfile&gt; &lt;inputfile&gt;</tt> should
 * produce a Postscript file describing an HMM.
 */
class HmmDrawerDot<H extends Hmm<?>> {
	public double minimumAij = 0.0;
	public double minimumPi = 0.0;
	protected NumberFormat probabilityFormat;

	/**
	 * This class converts an HMM to a dot file.
	 */
	public HmmDrawerDot() {
		probabilityFormat = NumberFormat.getInstance();
		probabilityFormat.setMaximumFractionDigits(2);
	}

	protected String convert(H hmm, double minimumAij, boolean simplify) {
		this.minimumAij = minimumAij;
		String s = beginning();

		if(simplify){
			s += simplifiedTransitions(hmm);
		}else{
			s += transitions(hmm);
		}
		s += states(hmm);

		return s + ending();
	}
	
	protected String convert(H hmm, double minimumAij) {
		this.minimumAij = minimumAij;
		String s = beginning();

		s += transitions(hmm);
		s += states(hmm);

		return s + ending();
	}
	protected String convertWODivision(H hmm, double minimumAij, List<ArrayList<Integer>> inState, List<ArrayList<Double>> inStateProbs) {
		this.minimumAij = minimumAij;
		String s = beginning();

		s += transitions(hmm);
		s += statesWODivision(hmm, inState, inStateProbs);

		return s + ending();
	}

	protected String convertSub(H hmm, double minimumAij, int depth) {
		this.minimumAij = minimumAij;
		String s = "";
		for (int i = 0; i < depth; i++){
			s += "\t";
		}
		s += beginning();

		s += transitionsSub(hmm, depth);
		s += statesSub(hmm, depth);
		
		for (int i = 0; i < depth; i++){
			s += "\t";
		}
		return s + ending();
	}


	protected String beginning() {
		return "digraph G {\n";
	}

	protected String transitions(Hmm<?> hmm) {
		String s = "";

		for (int i = 0; i < hmm.nbStates(); i++) {
			for (int j = 0; j < hmm.nbStates(); j++) {
				if (hmm.getAij(i, j) >= minimumAij && hmm.getAij(i, j) >= 0.000001)
					s += "\t" + i + " -> " + j + " [label=" + probabilityFormat.format(hmm.getAij(i, j)) + "];\n";
			}
		}
		return s;
	}
	protected String transitionsSub(Hmm<?> hmm, int depth) {
		String s = "";

		for (int i = 0; i < hmm.nbStates(); i++) {
			for (int j = 0; j < hmm.nbStates(); j++) {
				if (hmm.getAij(i, j) >= minimumAij && hmm.getAij(i, j) >= 0.001)
					for (int k = 0; k < depth; k++){
						s += "\t";
					}
					s += "\t" + i + " -> " + j + " [label=" + probabilityFormat.format(hmm.getAij(i, j)) + "];\n";
			}
		}
		return s;
	}
	
	protected String simplifiedTransitions(Hmm<?> hmm) {
		String s = "";
		
		/* find max output */
		for (int i = 0; i < hmm.nbStates() ; i++) { 
			double rowMax = 0;
			int colIndex = 0;
			for (int j = 0; j < hmm.nbStates(); j++) { 
				if(i==j && hmm.getAij(i, j) > minimumAij){
					s += "\t" + i + " -> " + j + " [label=" + probabilityFormat.format(hmm.getAij(i, j)) + "];\n";
				}else if (hmm.getAij(i, j) > rowMax){
					colIndex = j;
					rowMax = hmm.getAij(i, j);
				}
			}
			if(rowMax!=0){
				s += "\t" + i + " -> " + colIndex + " [label=" + probabilityFormat.format(rowMax) + "];\n";
//				System.out.println("\t" + i + " -> " + colIndex + " [label=" + probabilityFormat.format(rowMax) + "];");
			}
		}
		
		/* find max input */
		for (int j = 1; j < hmm.nbStates(); j++) { /* j != 0 */
			double colMax = 0;
			int rowIndex = 0;
			for (int i = 0; i < hmm.nbStates(); i++) {
				if (hmm.getAij(i, j) > colMax && i!=j) {
					rowIndex = i;
					colMax = hmm.getAij(i, j);
				}
			}
			if(!s.contains("\t" + rowIndex + " -> " + j) && 0!=colMax){
				s += "\t" + rowIndex + " -> " + j + " [label=" + probabilityFormat.format(colMax) + "];\n";
//				System.out.println("\t" + rowIndex + " -> " + j + " [label=" + probabilityFormat.format(colMax) + "];");
			}
		}
		
		return s;
	}


	protected String states(H hmm) {
		String s = "";

		for (int i = 0; i < hmm.nbStates(); i++) {
			
			s += "\t" + i + " [";

			if (hmm.getPi(i) >= minimumPi && hmm.getPi(i) >= 0.001) {
				s += "shape=rectangle, label=\"" + i + " - Pi= " + probabilityFormat.format(hmm.getPi(i)) + " - " + opdfLabel(hmm, i) + "\"";
			} else {
				s += "shape=rectangle, label=\"" + i + " - " + opdfLabel(hmm, i) + "\"";
			}

			s += "];\n";
		}

		return s;
	}
	protected String statesWODivision(H hmm, List<ArrayList<Integer>> inState, List<ArrayList<Double>> inStateProbs) {
		String s = "";

		for (int i = 0; i < hmm.nbStates(); i++) {
			ArrayList<Integer> activities = inState.get(i);
			
			s += "\t" + i + " [";

			if (hmm.getPi(i) >= minimumPi && hmm.getPi(i) >= 0.001) {
				s += "shape=rectangle, label=\"" + i + " - Pi= 1 - [ ";
				for (Integer activityIndex : activities) {
					s += SampleData.values()[activityIndex] + ",";
				}
				s += "]\"";
			} else {
				// s += "shape=rectangle, label=\"" + i + " - " + opdfLabelWODivision(hmm, i) + "\"";
				if (i == 2) {
					// end state
					s += "shape=rectangle, label=\"2 - [ end ";
				} else {
					s += "shape=rectangle, label=\"" + i + " - [ ";
					int index = 0;;
					for (Integer activityIndex : activities) {
						index = activities.indexOf(activityIndex);
						s += SampleData.values()[activityIndex] + "(" + probabilityFormat.format(inStateProbs.get(i).get(index)) + "),";
					}
				}
				s += "]\"";
			}

			s += "];\n";
		}

		return s;
	}
	protected String statesSub(H hmm, int depth) {
		String s = "";

		for (int i = 0; i < hmm.nbStates(); i++) {
			for (int j = 0; j < depth; j++){
				s += "\t";
			}
			s += "\t" + i + " [";

			if (hmm.getPi(i) >= minimumPi && hmm.getPi(i) >= 0.001) {
				s += "shape=rectangle, label=\"" + i + " - Pi= " + probabilityFormat.format(hmm.getPi(i)) + " - " + opdfLabel(hmm, i) + "\"";
			} else {
				s += "shape=rectangle, label=\"" + i + " - " + opdfLabel(hmm, i) + "\"";
			}

			s += "];\n";
		}

		return s;
	}

	protected String opdfLabel(H hmm, int stateNb)
	{
		Opdf<?> states = hmm.getOpdf(stateNb);		
		String str = "[ " + states.toString() + " ]";
		return str;
	}
	

	protected String ending() {
		return "}\n";
	}

	/**
	 * Writes a dot file depicting the given HMM.
	 *
	 * @param hmm
	 *            The HMM to depict.
	 * @param filename
	 *            The resulting 'dot' file filename.
	 */
	public void write(H hmm, String filename, double minimumAij, boolean simplify) throws IOException {
		FileWriter fw = new FileWriter(filename);
		fw.write(convert(hmm, minimumAij, simplify));
		
		fw.close();
	}
	public void writeSub(H hmm, String filename, double minimumAij, List<ArrayList<Integer>> inState, List<ArrayList<Double>> inStateProbs) throws IOException {
		FileWriter fw = new FileWriter(filename);
		fw.write(convertWODivision(hmm, minimumAij, inState, inStateProbs));
		
		fw.close();
	}

	public void print(H hmm, double minimumAij) {
		System.out.println(convert(hmm, minimumAij));
	}

	public void printSub(H hmm, double minimumAij, int depth) {
		// String tabString = "";
		// for (int i = 0; i < depth; i++) {
		// 	tabString += "\t";
		// }
		System.out.println(convertSub(hmm, minimumAij, depth));
	}
	
	public void write(H hmm, String filename) throws IOException {
		FileWriter fw = new FileWriter(filename);
		fw.write(convert(hmm, minimumAij));
		fw.close();
	}

	public void print(H hmm) {
		System.out.println(convert(hmm, minimumAij));
	}
}
