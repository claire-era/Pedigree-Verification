package preprocess;

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
	private static ArrayList<GenotypeTable> listOfGenos;
	private static GenotypeTable genos;

	public PreprocessHMP(String[] args) {
		PreprocessHMP.checkInputs(args);
	}

	private static String[][] readManualScores(String filename) throws IOException {
		File file = new File(filename);
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(file));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		String st;
		while ((st = br.readLine()) != null) {
			String s[] = st.split("\t");
//			System.out.println(s[0] + s[1]);
		}
		
		
		
		for(int taxa = 0; taxa <  genos.numberOfTaxa(); taxa++) {
			byte[] gene_arr = genos.genotypeAllSites(taxa); //returns genotype of taxa
//			Arrays.equals(gene_arr, gene_arr);
			
			System.out.print(genos.taxaName(taxa));
			System.out.print("\n");
			for(int j = 0 ;j < gene_arr.length; j++) {
				System.out.print(gene_arr[j]);
			}
			
			System.out.print("\n");
			
			
			
		}

		
		return null;
	}

	private static void checkInputs(String[] args) {
		try {
			if (args.length == 0) {
				System.out.println(
						"No arguments. \nPlease supply .hmp.txt file and pedigree file separated by space.\nFormat: java -jar PedVer.jar <input.hmp.txt> <pedigree_file.txt> <cut-off>");
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
		System.out.println("In " + hmp_filename + ":");
		return genos;
	}

	private static boolean init(GenotypeTable genos) throws IOException {
		boolean val = true;
		String[][] arr = PreprocessHMP.readManualScores("D:\\sp-workspace\\pevar-v0.1\\src\\input_file_4_preprocess.txt");
		// for every taxa in hapmap, check

		return val;
	}

	private static void start(String hmpFile, String pedFile, double cutOff) throws InterruptedException {
		try {
			// READ FILE
			genos = PreprocessHMP.ReadHMPFile(hmpFile);
			boolean value = PreprocessHMP.init(genos);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public static void main(String[] args) {
		PreprocessHMP p = new PreprocessHMP(args);
	}
	
}


