package pedigreeVerification;

import java.awt.FileDialog;
import java.awt.Frame;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

public class Step1 {
	private static ArrayList<GenotypeTable> listOfGenos;
	private static GenotypeTable genos;

	public Step1(String[] args) {
		Step1.checkInputs(args);
	}

	private static ArrayList<GenotypeTable> ReadFile(String hmpFile, String pedFile) { // Read file contents
		ArrayList<GenotypeTable> listOfGenos = PedigreeFileInfo.getFilteredGenotypeTable(hmpFile, pedFile);
		return listOfGenos;
	}

	private static GenotypeTable ReadHMPFile(String hmp_filename) {
		genos = ImportUtils.readFromHapmap(hmp_filename);
		return genos;
	}

	private static ArrayList<Integer> GetPolymorphicMarkers(GenotypeTable genos) {
		ArrayList<Integer> polyMarkers = new ArrayList<>();
		int hetCount, site = 0;
		int numOfSites = genos.numberOfSites(); // for 0-indexing
		int[] hcArr = new int[numOfSites];

		for (site = 0; site < numOfSites; site++) {
			hetCount = genos.heterozygousCount(site);
			hcArr[site] = hetCount;
		}

		// get indices of max count of polymorphism
		int max = hcArr[0];
		for (int i = 1; i < hcArr.length; i++)
			if (hcArr[i] > max)
				max = hcArr[i];

		// find all indices of polymorphic markers matching the max value
		for (site = 0; site < numOfSites; site++)
			if (genos.isPolymorphic(site) && max == genos.heterozygousCount(site))
				polyMarkers.add(site);

		return polyMarkers; // returns index of the location of the polymorphic marker
	}

	private static ArrayList<Integer> getInitSolution(ArrayList<Integer> polymarkers) {
		ArrayList<Integer> taxa = new ArrayList<Integer>();

		int psite = polymarkers.get(0); // get the max site
		for (int taxon = 0; taxon < genos.numberOfTaxa(); taxon++)
			if (genos.isHeterozygous(taxon, psite))
				taxa.add(taxon); // adds index of a sample included in the initial solution

		return taxa;
	}

	private static String createVector(double[] observed, double[] expected) {
		// CREATE VECTOR OF (OBSERVED,EXPECTED) values
		String res = "c(";
		for (int a = 0; a < observed.length; a++)
			res = res + Double.toString(observed[a]) + ",";

		for (int a = 0; a < expected.length; a++) {
			if (a == (expected.length - 1))
				res = res + Double.toString(expected[a]) + ")";
			else
				res = res + Double.toString(expected[a]) + ",";
		}
		// -----------------
		return res;
	}

	private static Double getPvalue(Rengine re, double[] observed, double[] expected) {
		String res = createVector(observed, expected);
		Double pval = null;
		if (!re.waitForR()) {
			System.out.println("Cannot load R");
			return null;
		}
		try {
			// ----- evaluate string
			String sample_m = "matrix(" + res + ",nrow = 2)";
			re.eval("observed=" + sample_m);
			re.eval("result=fisher.test(observed, alternative = \"greater\")");
			pval = re.eval("result$p.value").asDouble();
		} catch (Exception e) {
		}
		return pval;
	}

	private static ArrayList<Integer> GetCrossType(GenotypeTable genos) {
		ArrayList<Integer> ct = new ArrayList<Integer>(); // CROSS TYPE / ALLELE COMBINATION after function
		ArrayList<Integer> polymarkers = GetPolymorphicMarkers(genos);
		ArrayList<Integer> initTaxaSol = getInitSolution(polymarkers); // generated initial solution after good marker
																		// selection //taxa indices of initial solution
		Rengine re = new Rengine(null, false, null);
//		int psite;
		int nCountSites = genos.numberOfSites();
//		int pCountSites = polymarkers.size();
		int initSolTaxa = initTaxaSol.size();
		double p, q, X, Y, Z;
		int N;
		double[] expected1 = { 0, initSolTaxa, 0 }; // ASSUMPTION: ALL SAMPLES ARE HET (Segregation: AA x aa) 100% het
		double[] expected2 = { initSolTaxa / 2, initSolTaxa / 2, 0 }; // ASSUMPTION: HALF ARE HET AND HALF ARE HOMO
																		// DOMINANT (Segregation: AA x Aa) 50% self
		
		for (int site = 0; site < nCountSites; site++) {
			p = genos.majorAlleleFrequency(site);
			q = genos.minorAlleleFrequency(site);
			N = genos.numberOfTaxa();
			Y = genos.heterozygousCount(site); // HETEROZYGOUS
			X = ((2 * N * p) - Y) / 2; // HOMOZYGOUS DOM
			Z = ((2 * N * q) - Y) / 2; // HOMOZYGOUS REC

			double[] observed = { X, Y, Z };

			double pval_type1 = Step1.getPvalue(re, observed, expected1);
			double pval_type2 = Step1.getPvalue(re, observed, expected2);
			
			int type = (pval_type1 > pval_type2) ? 1 : 2;
			ct.add(type);
		}
		
		re.end(); //END R ENGINE INSTANCE AFTER PERFORMING FISHER'S EXACT TEST ON ALL SITES/SNPS
		return ct;
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

	public static void ProgramEnded() {
//		System.out.println("... Done.");
//		Thread.sleep(500);
		long startTime = System.nanoTime();
		long endTime = System.nanoTime();
		System.out.println("Took " + (endTime - startTime) + " ns");
	}

	private static void start(String hmpFile, String pedFile, double cutOff) throws InterruptedException {
		System.out.println("Cut-off: " + cutOff);

		try {
			// READ FILE
//			Step1.listOfGenos = ReadFile(hmpFile, pedFile);
			// GET KEYS FROM CORRESPONDING PEDIGREE FILE
//			ArrayList<String> keys = PedigreeFileInfo.getKeys();

			// IMPLEMENT STEP 1
			genos = Step1.ReadHMPFile(hmpFile);
			Step1.GetCrossType(genos);

			ProgramEnded();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
