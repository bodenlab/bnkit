package bn.example;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Iterator;
import java.util.Map;

import dat.Domain;
import dat.EnumVariable;
import dat.Enumerable;
import bn.prob.EnumDistrib;
import bn.alg.CGTable;
import bn.alg.Query;
import bn.alg.VarElim;
import bn.node.CPT;
import bn.BNet;
import bn.Predef;
/**
 * used for secondary structure prediction
 * with algiment 
 * @author wangyufei
 */

public class SecondaryStructurePredictor {
	
	final static private int windowSize = 11;
	final static private String fileName = "data/points.csv";
	final static private double trainDataRadio = 0.5f;
	static private List<String> testAmiodData;
	static private List<String> testSecondaryData;
	
	public static void main(String[] args) {
		System.out.println(NaiveBayes());
	}
	
	public static double NaiveBayes() {
		EnumVariable classNode = Predef.SecondaryStructure("secondary structure class");
		BNet bn = new BNet();
		
		//train
		List<EnumVariable> features = new LinkedList<EnumVariable>();
		initProcessData("data/result.out");
		List<Character[]> training = loadData("data/points.csv",windowSize);
		
		for(int i = 0; i < windowSize; i++) {
			features.add(Predef.AminoAcid("feature" + String.valueOf(i)));
		}
		double[][][] trainResult = learn(training,features,classNode);
		
		CPT head = new CPT(classNode);
		EnumDistrib classDistrib = new EnumDistrib(classNode.getDomain());
		
		double[] classCount = new double[classNode.size()];
		Arrays.fill(classCount, 0.0);
		
		Iterator<Character[]> iterator = training.iterator();
		Character[] record;
		while(iterator.hasNext()) {
			record = iterator.next();
			classCount[classNode.getIndex(record[windowSize])] ++;
		}
		
		classDistrib.set(classCount);
		head.put(classDistrib);
		
		bn.add(head);
		List<CPT> featureList = new LinkedList<CPT>();
		for(int i = 0; i < windowSize; i++) {
			featureList.add(new CPT(features.get(i),classNode));
			for(int j = 0; j < classNode.size(); j++){
				EnumDistrib distrib = new EnumDistrib(features.get(i).getDomain());
				distrib.set(trainResult[i][j]);
				featureList.get(i).put(distrib,classNode.getDomain().get(j));
			}
			bn.add(featureList.get(i));
		}
		
		VarElim ve = new VarElim();
        ve.instantiate(bn);
        Query q = null;
        double rightCase = 0;
        double totalCase = 0;
        Iterator<String> iter1 = testAmiodData.iterator();
        Iterator<String> iter2 = testSecondaryData.iterator();
        String amiodData;
        String secondaryData;
        EnumDistrib distrib;
        while(iter1.hasNext() && iter2.hasNext()) {
        	amiodData = iter1.next();
        	secondaryData = iter2.next();
        	for(int i = 0; i < amiodData.length() - windowSize + 1; i++) {
        		for(int j = 0; j < featureList.size(); j++){
        			featureList.get(j).setInstance((Character)amiodData.charAt(i + j));
        		}
        		q = ve.makeQuery(classNode);
    			CGTable r0 = (CGTable)ve.infer(q);
    			distrib = (EnumDistrib)r0.query(classNode);
                totalCase ++;
                if(classNode.getIndex((Character)secondaryData.charAt(i + windowSize / 2)) 
                		== distrib.getMaxIndex()) {
                	rightCase ++;
                }
			}
			
        }
		File toDelete = new File(fileName);
		if(toDelete.exists() && toDelete.isFile()) {
			toDelete.delete();
		}
		return rightCase / totalCase;
	}
	
	public static void initProcessData(String initDataFileName) {
		BufferedReader br = null;
		BufferedWriter bw = null;
		String amiod = "";
		String secondary = "";
		String[] data;
		List<String> initData = new LinkedList<String>();
		List<String> initTrainingData ;
		List<String> initTestData;
		try{
			br = new BufferedReader(new FileReader(initDataFileName));
			bw = new BufferedWriter(new FileWriter(new File(fileName))); 
			while ((amiod = br.readLine()) != null) { 
				secondary = br.readLine();
				initData.add(amiod + ";" + secondary);
			}
			Collections.shuffle(initData);
			int cutOffIndex = (int)(initData.size() * trainDataRadio);
			initTrainingData = initData.subList(0, cutOffIndex);
			initTestData = initData.subList(cutOffIndex, initData.size());
			
			Iterator<String> iter = initTrainingData.iterator();
			 
			while(iter.hasNext()) {
				data = iter.next().split(";");
				amiod = data[0];
				secondary = data[1];
				
				for(int i = 0; i < amiod.length() - windowSize + 1; i++) {
					bw.write(amiod.substring(i, i + windowSize) + "," + secondary.charAt(i + windowSize / 2));
					bw.newLine();  
	                bw.flush();  
				}
			}
			
			testAmiodData = new LinkedList<String>();
			testSecondaryData = new LinkedList<String>();
			
			iter = initTestData.iterator();
			while(iter.hasNext()) {
				data = iter.next().split(";");
				testAmiodData.add(data[0]);
				testSecondaryData.add(data[1]);
			}
			
		}catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
					bw.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public static double[][][] learn(List<Character[]> data, List<EnumVariable> features, EnumVariable parent) {
		double[][][] result = new double[features.size()][parent.size()][];
		for(int i = 0; i < features.size(); i++) {
			for(int j = 0; j < parent.size(); j++) {
				result[i][j] = new double[features.get(i).size()];
				for(int x = 0; x < features.get(i).size(); x++) {
					result[i][j][x] = 0;
				}
			}
		}
		
		Iterator<Character[]> iterator = data.iterator();
		while(iterator.hasNext()) {
			Character[] record = iterator.next();
			int parentIdx = parent.getIndex(record[features.size()]);
			for(int i = 0; i < features.size(); i++) {
				int valueIdx = features.get(i).getIndex(record[i]);
				result[i][parentIdx][valueIdx] ++;
			}
		}
		
		return result;
	}
	
	/**
	 * load data from scv file
	 * @param fileName the data file
	 * @param attributeSize the attribute number
	 * @return the list of all attribute and label for class
	 */
	public static List<Character[]> loadData(String fileName, int attributeSize) {
		List<Character[]> resultList = new LinkedList<Character[]>();
		String line = "";
		BufferedReader br = null;
		try{
			br = new BufferedReader(new FileReader(fileName));
		
			while ((line = br.readLine()) != null) { 
				// use comma as separator
				String[] content = line.split(",");
				Character[] data  = new Character[attributeSize + 1];
				for(int i = 0; i < attributeSize; i++) {
					data[i] = content[0].charAt(i);
				}
				data[attributeSize] = content[1].charAt(0);
				resultList.add(data);
			}
			
		}catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		return resultList;
	}
	
}
