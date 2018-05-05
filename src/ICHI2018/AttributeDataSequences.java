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

public class AttributeDataSequences {
	private int id;
	private int[] attributes = new int[18];

	
	public AttributeDataSequences(int id) {
		this.id = id;
	
	}

	public void setAttribute(int index, int value){
		this.attributes[index] = value;
	}

	public int getAttribute(int index) {
		return this.attributes[index];
	}
	public int getId() {
		return this.id;
	}
}