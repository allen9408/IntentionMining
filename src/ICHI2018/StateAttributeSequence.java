package ICHI2018;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import HMM.Observation.ObservationDiscrete;

public class StateAttributeSequence {
	private int attributeNumber = 18;
	private ArrayList<ArrayList<Integer>> attributeSequence = new ArrayList<ArrayList<Integer>>(attributeNumber);
	private int StateIndex;

	public StateAttributeSequence(int index) {
		StateIndex = index;
		for (int i =0; i < attributeNumber; i++) {
			ArrayList<Integer> list = new ArrayList<Integer>();
			attributeSequence.add(list);
		}
	}

	public void addAttributeContent(int attributeIndex, int attributeValue) {
		ArrayList<Integer> attributesContent = attributeSequence.get(attributeIndex);
		attributesContent.add(attributeValue);
	}

	public int getAttributeValueNumber(int attributeIndex) {
		ArrayList<Integer> attributesContent = attributeSequence.get(attributeIndex);
		ArrayList<Integer> tmp = new ArrayList<Integer>();
		int number = 0;
		for (int value : attributesContent) {
			if (!tmp.contains(value)) {
				number++;
				tmp.add(value);
			}
		}
		return number;
	}

	public ArrayList<Integer> getAttributeValueSequence(int attributeIndex) {
		ArrayList<Integer> attributesContent = attributeSequence.get(attributeIndex);
		ArrayList<Integer> tmp = new ArrayList<Integer>();
		for (int value : attributesContent) {
			if (!tmp.contains(value)) {
				tmp.add(value);
			}
		}
		return tmp;
	}

	public double getAttributeProbability(int attributeIndex, int value) {
		ArrayList<Integer> attributesContent = attributeSequence.get(attributeIndex);
		Collections.sort(attributesContent);
		int appear = attributesContent.lastIndexOf(value) - attributesContent.indexOf(value) + 1;
		double prob = appear;
		return prob;
	}
	
	public double size(int attributeIndex) {
		ArrayList<Integer> attributesContent = attributeSequence.get(attributeIndex);
		return attributesContent.size();
	}

}