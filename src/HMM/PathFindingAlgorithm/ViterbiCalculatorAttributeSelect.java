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

package HMM.PathFindingAlgorithm;

import java.util.Iterator;
import java.util.List;
import java.util.Hashtable;

import HMM.DataType.Hmm;
import HMM.Observation.Observation;
import ICHI2018.AttributeDataSequences;
import ICHI2018.StateAttributeSequence;


/**
 * This class can be used to compute the most probable state sequence matching
 * a given observation sequence (given an HMM).
 */
public class ViterbiCalculatorAttributeSelect
{	
	/*
	 * The psy and delta values, as described in Rabiner and Juand classical
	 * papers.
	 */
	private double[][] delta; 
	private int[][] psy;
	private int[] stateSequence;
	private double lnProbability;
	
	
	/**
	 * Computes the most likely state sequence matching an observation
	 * sequence given an HMM.
	 *
	 * @param hmm A Hidden Markov Model;
	 * @param oseq An observations sequence.
	 */
	public <O extends Observation> 
	ViterbiCalculatorAttributeSelect(List<? extends O> oseq, Hmm<O> hmm, Hashtable<Integer, AttributeDataSequences> attributeSet, StateAttributeSequence[][] transitionAttributeSet, int idNum, List<StateAttributeSequence> stateAttributeSet,
		List<Integer> ignoredAttributes, int currentIgnoredAttribute)
	{
		if (oseq.isEmpty())
			throw new IllegalArgumentException("Invalid empty sequence");
		
		delta = new double[oseq.size()][hmm.nbStates()];
		psy = new int[oseq.size()][hmm.nbStates()];
		stateSequence = new int[oseq.size()];
		
		for (int i = 0; i < hmm.nbStates(); i++) {
			delta[0][i] = -Math.log(hmm.getPi(i)) - 
			Math.log(hmm.getOpdf(i).probability(oseq.get(0)));
			psy[0][i] = 0;
		}
		
		Iterator<? extends O> oseqIterator = oseq.iterator();
		if (oseqIterator.hasNext())
			oseqIterator.next();
		
		int t = 1;
		// Get attribute Sequence
		AttributeDataSequences attributeSequence = attributeSet.get(idNum);
		while (oseqIterator.hasNext()) {
			O observation = oseqIterator.next();
			
			for (int i = 0; i < hmm.nbStates(); i++)
				computeStep(hmm, observation, t, i, attributeSequence, transitionAttributeSet, stateAttributeSet, ignoredAttributes, currentIgnoredAttribute);
			
			t++;
		}
		
		lnProbability = Double.MAX_VALUE;
		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisProbability = delta[oseq.size()-1][i];
			
			if (lnProbability > thisProbability) {
				lnProbability = thisProbability;
				stateSequence[oseq.size() - 1] = i;
			}
		}
		lnProbability = -lnProbability;
		
		for (int t2 = oseq.size() - 2; t2 >= 0; t2--)
			stateSequence[t2] = psy[t2+1][stateSequence[t2+1]];
	}
	
	
	/*
	 * Computes delta and psy[t][j] (t > 0) 
	 */
	private <O extends Observation> void
	computeStep(Hmm<O> hmm, O o, int t, int j, AttributeDataSequences attributeSequence, StateAttributeSequence[][] transitionAttributeSet, List<StateAttributeSequence> stateAttributeSet,
		List<Integer> ignoredAttributes, int currentIgnoredAttribute) 
	{
		double minDelta = Double.MAX_VALUE;
		int min_psy = 0;
		int[] typeNum = {4, 2, 3, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
		
		for (int i = 0; i < hmm.nbStates(); i++) {
			StateAttributeSequence tmpSequence = transitionAttributeSet[i][j];
			// smmothing
			double[] attrProbs = new double[18];
			if (transitionAttributeSet[i][j].size(0) == 0) {
				// no attribute record from transition i to j
				for (int k = 0; k < 18; k++) {
					attrProbs[k] = 1.0 / typeNum[k];
				}
			} else {
				// Already has attribute record, smoothing for each attribute value
				for (int k = 0; k < 18; k++) {
					// get ID's attribute
					int attributeValue = attributeSequence.getAttribute(k);
					int valueNumber = transitionAttributeSet[i][j].getAttributeValueNumber(k);
					// smoothing
					if (transitionAttributeSet[i][j].getAttributeProbability(k, attributeValue) == 0) {
						attrProbs[k] = 0.05/typeNum[k];
					} else {
						double prob = transitionAttributeSet[i][j].getAttributeProbability(k, attributeValue)/transitionAttributeSet[i][j].size(k);
						attrProbs[k] = prob * (1 - 0.05) + 0.05/typeNum[k];
					}
				}

			}

			double attributeLnProb = 0.0;
			for (int k = 0; k < 18; k++) {
				if (ignoredAttributes.contains(k) || currentIgnoredAttribute == k) {
					continue;
				}
				attributeLnProb += Math.log(attrProbs[k]);
			}
			double thisDelta = delta[t-1][i] - Math.log(hmm.getAij(i, j)) - attributeLnProb/18;
			
			if (minDelta > thisDelta) {
				minDelta = thisDelta;
				min_psy = i;
			}
		}

		// Add emmision attribute probability
		StateAttributeSequence emmisionAttributeSequence = stateAttributeSet.get(j);
		double[] emmisionProbs = new double[18];
		if (emmisionAttributeSequence.size(0) == 0) {
			for (int k = 0; k < 18; k++) {
				emmisionProbs[k] = 1.0 / typeNum[k];
			}
		} else {
			for (int k = 0; k < 18; k++) {
				int attributeValue = attributeSequence.getAttribute(k);
				// smoothing
				if (emmisionAttributeSequence.getAttributeProbability(k, attributeValue) == 0) {
					emmisionProbs[k] = 0.05/typeNum[k];
				} else {
					double prob = emmisionAttributeSequence.getAttributeProbability(k, attributeValue)/emmisionAttributeSequence.size(k);
					emmisionProbs[k] = prob * (1 - 0.05) + 0.05/typeNum[k];					
				}
			}
		}
		double emmisionAttributeProb = 0.0;
		for (int k = 0; k < 18; k++) {
			if (ignoredAttributes.contains(k) || currentIgnoredAttribute == k) {
				continue;
			}
			emmisionAttributeProb += Math.log(emmisionProbs[k]);
		}
		
		delta[t][j] = minDelta - Math.log(hmm.getOpdf(j).probability(o)) - emmisionAttributeProb/18;
		psy[t][j] = min_psy;
	}
	
	
	/**
	 * Returns the neperian logarithm of the probability of the given
	 * observation sequence on the most likely state sequence of the given
	 * HMM.
	 *
	 * @return <code>ln(P[O,S|H])</code> where <code>O</code> is the given
	 *         observation sequence, <code>H</code> the given HMM and 
	 *         <code>S</code> the most likely state sequence of this observation
	 *         sequence given this HMM.
	 */
	public double lnProbability()
	{
		return lnProbability;
	}
	
	
	/**
	 * Returns a (clone of) the array containing the computed most likely
	 * state sequence.
	 *
	 * @return The state sequence; the i-th value of the array is the index
	 *         of the i-th state of the state sequence.
	 */
	public int[] stateSequence() 
	{
		return stateSequence.clone();
	}
}
