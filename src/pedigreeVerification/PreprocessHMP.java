package pedigreeVerification;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import pedigreeVerification.Step1;

public class PreprocessHMP {
//	private static ArrayList<GenotypeTable> listOfGenos;
	private static GenotypeTable genos;
	private static boolean value;
//	private static boolean isF1 = false;
	private static String fname;

	private static ArrayList<Integer> scores_manual;

	public PreprocessHMP(String[] args) {
		PreprocessHMP.checkInputs(args);
	}

	private static boolean readManualScores(String filename, String hmpFile) throws IOException {
		File file = new File(filename);
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(file));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		String st;
		ArrayList<Integer> score = new ArrayList<Integer>();
		String subs = hmpFile.split(".hmp.txt")[0]; // get group name// group file name
		PreprocessHMP.fname = hmpFile;
		while ((st = br.readLine()) != null) {
			String s[] = st.split("\t");
//			System.out.println(s[0] + " "  +  s[1]+ " " + s[2]);
//			System.out.println(s[0].contains(subs));
			
//			System.out.println(s[0] + " " + s[2]);
			
//			if ((s[2].equals("F1"))) { // if not F1, most likely parent sya
//				PreprocessHMP.setF1(true);
//			}

			if (s[0].contains(subs)) {// compare to current groups
//				System.out.println(s[0] + " " + s[2]);
				score.add(Integer.parseInt(s[1]));
			}
		}
		// CLASSIFICATION OF STRINGS
		ArrayList<String> genotypeAllTaxa = new ArrayList<String>();
		for (int taxon = 0; taxon < genos.numberOfTaxa(); taxon++) {
			genotypeAllTaxa.add(genos.genotypeAsStringRow(taxon));
//			System.out.println(score.get(taxon) + " " + genotypeAllTaxa.get(taxon));
		}

		// get classifications.
		int class_i = 0;
		ArrayList<String> classif = new ArrayList<String>(); // contains unique values
		ArrayList<Integer> index_class = new ArrayList<Integer>(); // contains classifcation assignment values
		for (int taxon = 0; taxon < genotypeAllTaxa.size(); taxon++) {
			String seq = genotypeAllTaxa.get(taxon);
			if (!classif.contains(seq)) { // if value is inside the list, check its index.?
				// assign new classification
				classif.add(seq);
				index_class.add(class_i);
				class_i++;
			} else { // else if seq already has classification
				int curr_index_class = classif.indexOf(seq);
				index_class.add(curr_index_class);
			}
		}

		// for every new classification encountered, set reference score.
		ArrayList<Integer> ref = new ArrayList<Integer>();
		for (int c = 0; c < classif.size(); c++) {
			for (int t = 0; t < genos.numberOfTaxa(); t++) {
				if (index_class.get(t) == c) {
					ref.add(score.get(t));
					break;
				}
			}
		}

//		for(int i =0; i < ref.size(); i++) {
//			System.out.println(ref.get(i));
//		}
		int curr_score = 0;
		int curr_class = 0;
		boolean equal = false;
		for (int j = 0; j < genos.numberOfTaxa(); j++) { // MAP CLASSIFICATION TO INDEX OF C
			equal = false;
			curr_score = score.get(j);
			curr_class = index_class.get(j);
			if (ref.get(curr_class) == curr_score) { // if reference == current
				equal = true;
			} else
				equal = false;

			if (equal) {
				// meaning all taxa have equal classifications!!!
//				System.out.println("equal sis");
			} else {
//				System.out.println("may error sa data natin sis");
				break;
			}
		}

//		scores_manual = scores;
//		for(int i = 0; i < score.size(); i++) {
//			Step1.scores_manual[i] = score.get(i);
//		}

		PreprocessHMP.scores_manual = score;

		return equal;
	}

	private static void checkInputs(String[] args) {
		try {
			if (args.length == 0) {
				System.out.println("No arguments. \nPlease supply .hmp.txt file and pedigree file separated by space.\nFormat: java -jar PedVer.jar <input.hmp.txt> <pedigree_file.txt> <cut-off>");
			} else if (args.length == 3) { // then there must be a provided cut-off
				start(args[0], args[1], Double.parseDouble(args[2]));
			} else if (args.length == 2) { // then default, get minimum
				start(args[0], args[1], 0.0);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static GenotypeTable ReadHMPFile(String hmp_filename) {
		genos = ImportUtils.readFromHapmap(hmp_filename);
//		System.out.println("In " + hmp_filename + ":");
		return genos;
	}

	private static boolean init(GenotypeTable genos, String hmpFile) throws IOException {
		boolean val = false;
		val = PreprocessHMP.readManualScores("D:\\sp-workspace\\pevar-v0.1\\src\\input_preprocess.txt", hmpFile);
		// for every taxa in hapmap, check

		return val;
	}

	private static void start(String hmpFile, String pedFile, double cutOff) throws InterruptedException {
		try {
			// READ FILE
			if((hmpFile.contains(".hmp.txt"))) {
				genos = PreprocessHMP.ReadHMPFile(hmpFile);
				boolean value = PreprocessHMP.init(genos, hmpFile);
				PreprocessHMP.setVal(value);
			}else {
				System.out.println("Invalid data input.\n");
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void setVal(boolean value) {
		PreprocessHMP.value = value;
	}
//
//	private static void setF1(boolean b) {
//		PreprocessHMP.isF1 = b;
//	}
	
	public String getFilename() {
		return PreprocessHMP.fname;
	}
	public boolean getVal() {
		return PreprocessHMP.value;
	}

	public ArrayList<Integer> getManualScores() {
		return PreprocessHMP.scores_manual;
	}

//	public boolean isF1() {
//		if ((PreprocessHMP.isF1)) {
//			return PreprocessHMP.isF1;
//		}
//		else
//			return false;
//	}

}
