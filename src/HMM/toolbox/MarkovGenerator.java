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

package HMM.toolbox;

import java.util.ArrayList;
import java.util.List;

import HMM.DataType.Hmm;
import HMM.Observation.Observation;

/**
 * Generates sequences of markovian observations given a HMM.
 */

public class MarkovGenerator<O extends Observation> {
	private final Hmm<O> hmm;
	private int stateNb;
	private int endNb;
	private double seqLogLikelihood = 0.0;
	private List<Integer> states = new ArrayList<Integer>();

	/**
	 * Initializes a Markovian generator.
	 *
	 * @param hmm
	 *            An Hidden Markov Model that perfectly models the sequences
	 *            generated by this object.
	 */
	public MarkovGenerator(Hmm<O> hmm) {
		if (hmm == null)
			throw new IllegalArgumentException("Invalid null HMM");

		this.hmm = hmm;
		newSequence();
	}

	/**
	 * Generates a new (pseudo) random observation.
	 *
	 * @return The generated observation.
	 */
	public O observation() {
		O o = hmm.getOpdf(stateNb).generate();
		double rand = Math.random();

		
		for (int j = 0; j < hmm.nbStates() - 1; j++)
			if ((rand -= hmm.getAij(stateNb, j)) < 0) {
				stateNb = j;
				states.add(j);
				if (hmm.getAij(stateNb, j) != 0)
					seqLogLikelihood += Math.log(hmm.getAij(stateNb, j));
				// System.out.print(stateNb + ", ");
				return o;
			}

		stateNb = hmm.nbStates() - 1;
		return o;
	}

	public O observation_scaled(int curLength, int avgLength) {
		O o = hmm.getOpdf(stateNb).generate();
		double rand = Math.random();

		double[] aij = new double[hmm.nbStates()];
		for (int i = 0; i < hmm.nbStates(); i++) {
			aij[i] = hmm.getAij(stateNb, i);
		}
		// scale for avg Length
		double para = 1.0;
		if (curLength < avgLength) {
			double probToEnd = aij[2];
			for (int i = 0; i < hmm.nbStates(); i++) {
				aij[i] += para * probToEnd * (1 - (curLength/avgLength)*curLength/avgLength) / (hmm.nbStates() - 1);
			}
			aij[2] = probToEnd * curLength/avgLength*curLength/avgLength / para;
		} else {
			for (int i = 0; i < hmm.nbStates(); i++) {
				if (i == 2) continue;
				aij[i] *= (1 - 5 * (curLength -  avgLength)/curLength/hmm.nbStates());
			}
			aij[2] += 5 * (1 - aij[2]) * ((curLength - avgLength)/curLength);
		}
		
		for (int j = 0; j < hmm.nbStates() - 1; j++)
			if ((rand -= aij[j]) < 0) {
				stateNb = j;
				states.add(j);
				if (hmm.getAij(stateNb, j) != 0)
					seqLogLikelihood += Math.log(hmm.getAij(stateNb, j));
				// System.out.print(stateNb + ", ");
				return o;
			}

		stateNb = hmm.nbStates() - 1;
		return o;
	}

	/**
	 * Generates a new (pseudo) random observation sequence and start a new one.
	 * 
	 * @param length
	 *            The length of the sequence.
	 * @return An observation sequence.
	 */
	public List<O> observationSequence(int length) {
		if (length <= 0)
			throw new IllegalArgumentException("Positive length required");

		ArrayList<O> sequence = new ArrayList<O>();
		while (length-- > 0)
			sequence.add(observation());
		newSequence();
		// System.out.println(sequence);
		return sequence;
	}

	/* random length */
	public List<O> observationSequence() {

		ArrayList<O> sequence = new ArrayList<O>();
		// do {
		// O o = observation();
		// sequence.add(o);
		// } while (!sequence.get(sequence.size() - 1).equals("end"));
		// while (stateNb != endNb){
		// O o = observation();
		// sequence.add(o);
		// }
		// System.out.print("State sequence: ");
		while (true) {
			
			O o = observation();
			sequence.add(o);
			// System.out.println(o);
			if (o.toString().equals("end")) {
				break;
			}
		}
		// System.out.println("");
		// sequence.add(observation());
		newSequence();
		return sequence;
	}
	/* random length modified for length*/
	public List<O> observationSequence_scaled(int avgLength) {

		ArrayList<O> sequence = new ArrayList<O>();
		int curLength = 0;
		while (true) {
			
			O o = observation_scaled(curLength, avgLength);
			sequence.add(o);
			// System.out.println(o);
			if (o.toString().equals("end")) {
				break;
			}
			curLength++;
		}
		// System.out.println("");
		// sequence.add(observation());
		newSequence();
		return sequence;
	}

	/**
	 * Finds a new state according to the initial (pi) probabilities of each
	 * state.
	 */
	public void newSequence() {
		double rand = Math.random(), current = 0.;

		for (int i = 0; i < hmm.nbStates() - 1; i++) {
			current += hmm.getPi(i);

			if (current > rand) {
				stateNb = i;
				return;
			}
		}

		stateNb = hmm.nbStates() - 1;

		current = 0;
		for (int i = 0; i < hmm.nbStates() - 1; i++) {
			current += hmm.getEpi(i);
			if (current > rand) {
				stateNb = i;
				return;
			}
		}
		endNb = hmm.nbStates() - 1;
	}

	/**
	 * Returns the state number of the current state.
	 *
	 * @return A state number.
	 */
	public int stateNb() {
		return stateNb;
	}

	public double logLikely() {
		return seqLogLikelihood;
	}

	public List<Integer> getStates() {
		return states;
	}
}
