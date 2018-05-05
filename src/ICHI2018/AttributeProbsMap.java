package ICHI2018;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Enumeration;
public class AttributeProbsMap {
	private List<Hashtable<Integer, Double>> attrprobsSet = new ArrayList<Hashtable<Integer, Double>>(18);
	
	public AttributeProbsMap() {
		for (int i = 0; i < 18; i ++) {
			attrprobsSet.add(new Hashtable<Integer, Double>());
		}
	}

	public void addProbToSet(int attributeIndex, int value, double prob) {
		Hashtable<Integer, Double> attrproSeq = attrprobsSet.get(attributeIndex);
		if (attrproSeq.containsKey(value)) {
			double oldProb = attrproSeq.get(value);
			double newProb = oldProb + prob;
			// System.out.println("New Prob: " + newProb);
			attrproSeq.replace(value, newProb);
		} else {
			attrproSeq.put(value, prob);
		}
		attrprobsSet.set(attributeIndex, attrproSeq);
	}

	public Hashtable<Integer, Double> getAttributeProbs(int attributeIndex) {
		return attrprobsSet.get(attributeIndex);
	}

	public void mergeAttribute(int index, Hashtable<Integer, Double> attr) {
		Hashtable<Integer, Double> attrproSeq = attrprobsSet.get(index);
		Enumeration<Integer> e = attr.keys();
		while(e.hasMoreElements()) {
			int value = e.nextElement();
			double prob = attr.get(value);
			if (attrproSeq.containsKey(value)) {
				double oldProb = attrproSeq.get(value);
				double newProb = oldProb + prob;
				// System.out.println("New Prob: " + newProb);
				attrproSeq.replace(value, newProb);
			} else {
				attrproSeq.put(value, prob);
			}
		}
		attrprobsSet.set(index, attrproSeq);
	}

	public List<Integer> getMostLikelyAttributes() {
		List<Integer> result = new ArrayList<Integer>();
		int index = 0;
		for (Hashtable<Integer, Double> attr : attrprobsSet) {
			// System.out.println("Attribute: " + index);
			int maxValue = 0;
			double maxProb = 0.0;
			Enumeration<Integer> e = attr.keys();
			while(e.hasMoreElements()) {
				int value = e.nextElement();
				double prob = attr.get(value);
				if (prob > maxProb) {
					maxValue = value;
					maxProb = prob;
				}
				// System.out.println(value + ", " + attr.get(value));
			}
			result.add(maxValue);
			index++;
		}
		return result;
	}

	public List<Integer> getEstimatedAttributes() {
		List<Integer> result = new ArrayList<Integer>();
		for (Hashtable<Integer,Double> attr: attrprobsSet) {
			double sum = 0;
			Enumeration<Integer> e = attr.keys();
			while(e.hasMoreElements()) {
				int value = e.nextElement();
				sum += attr.get(value);
			}
			double rand = Math.random() * sum;
			Enumeration<Integer> f = attr.keys();
			while(f.hasMoreElements()) {
				int attributeValue = f.nextElement();
				if ((rand -= attr.get(attributeValue)) < 0) {
					result.add(attributeValue);
					break;
				}
			}

		}
		return result;
	}
}