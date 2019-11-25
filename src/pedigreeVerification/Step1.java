package pedigreeVerification;

import java.util.ArrayList;

import org.rosuda.JRI.Rengine;
import java.lang.Math;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

public class Step1 {
	private static ArrayList<GenotypeTable> listOfGenos;
	private static GenotypeTable genos;
	private static int scores_auto[];
	private static int scores_auto_selfing[];
//	public static Integer[] scores_manual;

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
		boolean polyexistence = false;
		
		//IMPORTANT: 1 PASS THRU ALL SITES PER TAXA AND SEE IF THERE IS ANY POLYMORPHISM. ELSE RETURN NULL, AND CONCLUDE THAT FOR ALL SITES, SAMPLES MUST BE SELFING
		for (int s = 0; s < numOfSites; s++) {
			if(genos.isPolymorphic(s)) {
				polyexistence = true;
				break;
			}
		}
		
		if(polyexistence == false) return null;
		
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
			if (genos.isPolymorphic(site) && max == genos.heterozygousCount(site)) {
				polyMarkers.add(site);
			}

		return polyMarkers; // returns index of the location of the polymorphic marker

	}

	private static ArrayList<Integer> getInitSolution(ArrayList<Integer> polymarkers) {
		ArrayList<Integer> taxa = new ArrayList<Integer>();

		int psite = polymarkers.get(0); // get the max site
//		System.out.println(psite + " psite to ");
		for (int taxon = 0; taxon < genos.numberOfTaxa(); taxon++)
			if (genos.isHeterozygous(taxon, psite)) {
//				System.out.println(taxon);
				taxa.add(taxon); // adds index of a sample included in the initial solution
			}

		return taxa;
	}

	private static double sum(double[] exp) {
		double s = 0;
		for (int i = 0; i < exp.length; i++)
			s += exp[i];
		return s;
	}

	private static double[] prod(double[] exp, double[] obs) {
		// multiply (sum(obs) * exp ) / sum(exp) [vector mult]
//		double s = sum(obs) * exp;
		double s = sum(obs);
		double e = sum(exp);
//		double res = 0;
//		double arr[] = new double[obs.length];
		double res[] = new double[obs.length];
		for (int i = 0; i < obs.length; i++) {
			exp[i] = s * exp[i];
			res[i] = exp[i] / e;
		}
		// get ceiling of all expected values
		for (int i = 0; i < obs.length; i++) {
			exp[i] = Math.ceil(res[i]);
		}

		// if for all exp, 0 lahat...
		if (exp[0] == 0 && exp[1] == 0 && exp[2] == 0) {
			exp[0] = 1;
			exp[1] = 1;
			exp[2] = 1;
		}

		return exp;
	}

	private static String createVector(double[] observed, double[] expected) {

		expected = Step1.prod(expected, observed);
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
		String vector_of_obs_exp = createVector(observed, expected);
		Double pval = null;
		if (!re.waitForR()) {
			System.out.println("Cannot load R");
			return null;
		}
		try {
			// ----- evaluate string
			String sample_m = "matrix(" + vector_of_obs_exp + ",nrow = 3)"; // for each element inside the
																			// classification array. expected[0]..[2].
			re.eval("matrixof_exp_obs=" + sample_m);// store as an R object
			re.eval("result=fisher.test(matrixof_exp_obs)");// call fisher.test and store in result
			pval = re.eval("result$p.value").asDouble(); // get p-value from list of results
//			System.out.println(pval);
		} catch (Exception e) {
		}
//		System.out.print("\n");
		return pval;
	}

	private static int getMax(double[] arr) { // returns index of max value in an array
		double max = arr[0];
		int i = 0;
		int max_i = 0;
		for (i = 1; i < arr.length - 1; i++)
			if (arr[i] > max) {
				max = arr[i];
				max_i = i;
			}
//		System.out.println(i);
		return max_i;

	}

	private static double countAlleles(ArrayList<Integer> listOftaxa, int site, GenotypeTable genos, int type) {
		double majCount = 0;
		double minCount = 0;
		double hetCount = 0;
		byte majorAlleleInSite = genos.majorAllele(site);
		byte minorAlleleInSite = genos.minorAllele(site);

		for (int i = 0; i < listOftaxa.size(); i++) {
			int t = listOftaxa.get(i);
			if (genos.isHeterozygous(t, site))
				hetCount++;
			else if (majorAlleleInSite == genos.genotypeArray(t, site)[1])
				majCount++;
			else if (minorAlleleInSite == genos.genotypeArray(t, site)[1])
				minCount++;
		}
		if (type == 0)
			return majCount;
		else if (type == 1)
			return minCount;
		else
			return hetCount;
	}

	private static ArrayList<Integer> GetCrossType(GenotypeTable genos) {
		ArrayList<Integer> ct = new ArrayList<Integer>(); // CROSS TYPE / ALLELE COMBINATION after function
		ArrayList<Integer> polymarkers = GetPolymorphicMarkers(genos);
		if(polymarkers == null) return null;
		
		ArrayList<Integer> initTaxaSol = getInitSolution(polymarkers); // generated initial solution after good marker
																		// selection //taxa indices of initial solution

		Rengine re = new Rengine(null, false, null);
		int nCountSites = genos.numberOfSites();
		double initSolTaxa = initTaxaSol.size();
		double X, Y, Z;
	
		if(initSolTaxa == 0) return null;
		
		double[] expected1 = { initSolTaxa, 0, 0 };
		double[] expected2 = { 0, 0, initSolTaxa };
		double[] expected3 = { 0, initSolTaxa, 0 }; // ASSUMPTION: ALL SAMPLES ARE HET (Segregation: AA x aa) 100% het
		double[] expected4 = { initSolTaxa / 2, initSolTaxa / 2, 0 }; // ASSUMPTION: HALF ARE HET AND HALF ARE HOMO //
																		// DOMINANT (Segregation: AA x Aa) 50% self
		double[] expected5 = { 0, initSolTaxa / 2, initSolTaxa / 2 };
		double[] expected6 = { initSolTaxa / 4, initSolTaxa / 2, initSolTaxa / 4 };

//		for(int i = 0; i < polymarkers.size(); i++) {
//			System.out.println(polymarkers.get(i));
//		}

//		for(int i = 0; i < initTaxaSol.size(); i++) {
//			System.out.println(initTaxaSol.get(i));
//		}

		for (int site = 0; site < nCountSites; site++) {
			X = countAlleles(initTaxaSol, site, genos, 0); // 0 for major allele
			Z = countAlleles(initTaxaSol, site, genos, 1); // 1 for minor allele
			Y = countAlleles(initTaxaSol, site, genos, 2); // HETEROZYGOUS COUNT

			double[] observed = { X, Y, Z };

//			System.out.println(observed[0] + " " + observed[1] + " " + observed[2]);
//			System.out.println(expected1[0] + " " + expected1[1] + " " + expected1[2]);

			double p_t1 = Step1.getPvalue(re, observed, expected1);
			double p_t2 = Step1.getPvalue(re, observed, expected2);
			double p_t3 = Step1.getPvalue(re, observed, expected3);
			double p_t4 = Step1.getPvalue(re, observed, expected4);
			double p_t5 = Step1.getPvalue(re, observed, expected5);
			double p_t6 = Step1.getPvalue(re, observed, expected6);

			double[] p_arr = { p_t1, p_t2, p_t3, p_t4, p_t5, p_t6 };

			// classify: get max between 6 numbers
			int index_max = getMax(p_arr);
			index_max = index_max + 1; // for fixing 0-indexing, classification

			ct.add(index_max);
		}

		re.end(); // END R ENGINE INSTANCE AFTER PERFORMING fisher.test() ON ALL SITES/SNPS
		return ct;
	}

	private static byte[][][] InferGenotype(GenotypeTable genos, ArrayList<Integer> crossType) {
		// RETURNS BYTE 2D ARRAY OF F1 AND SELF EXPECTED CROSS
		byte[][][] exp = new byte[2][crossType.size()][2]; // return value
		byte majAllele;
		byte minAllele;
		for (int site = 0; site < crossType.size(); site++) { // for each SNP, with corresponding crosstype
			byte[] exp_f1 = new byte[2]; // exp0 holds higher/major bit,1 holds lower/min/het
			byte[] exp_self = new byte[2];
			int cross = crossType.get(site); // get crosstype per site
			majAllele = genos.majorAllele(site);
			minAllele = genos.minorAllele(site);
			// INCUR DOMINANCE / majority / minority
			if (cross == 1) { // MAJ X MAJ == 100% homo dom expected children
				exp_f1[0] = majAllele;
				exp_f1[1] = majAllele;

				exp_self[0] = majAllele;
				exp_self[1] = majAllele;
			} else if (cross == 2) { // MIN X MIN == 100% homo rec expected children
				exp_f1[0] = minAllele;
				exp_f1[1] = minAllele;

				exp_self[0] = minAllele;
				exp_self[1] = minAllele;
			} else if (cross == 3) { // MAJ X MIN == HET == 100% het exp
				exp_f1[0] = majAllele;
				exp_f1[1] = minAllele;

				exp_self[0] = majAllele;
				exp_self[1] = majAllele;
			} else if (cross == 4) { // MAJOR x HET == 50% het, 50% homo dom
				exp_f1[0] = majAllele; // doblehin mo nalang sis
				exp_f1[1] = GenotypeTableUtils.getDiploidValue(majAllele, minAllele);
//				System.out.println(NucleotideAlignmentConstants.getNucleotideIUPAC(exp_f1[1]));

				exp_self[0] = GenotypeTableUtils.getDiploidValue(majAllele, minAllele);
				;
				exp_self[1] = exp_self[0];
			} else if (cross == 5) { // HET X MINOR
				exp_f1[0] = GenotypeTableUtils.getDiploidValue(majAllele, minAllele);
				exp_f1[1] = minAllele;

				exp_self[0] = majAllele;
				exp_self[1] = majAllele;
			} else if (cross == 6) { // HET X HET
				exp_f1[0] = GenotypeTableUtils.getDiploidValue(majAllele, minAllele);
				exp_f1[1] = exp_f1[0];

				exp_self[0] = GenotypeTableUtils.getDiploidValue(majAllele, minAllele);
				;
				exp_self[1] = exp_self[0];
			}
			// STORE INSIDE RETURN ARRAY impose 0-indexing
			exp[0][site] = exp_f1;
			exp[1][site] = exp_self;
		}
		return exp;

	}

	private static int[] testSolution(GenotypeTable genos, byte[][] exp_f1, ArrayList<Integer> crossType) {
		int scores_auto[] = new int[genos.numberOfTaxa()];
		int count_1 = 0, count_0 = 0;
		for (int t = 0; t < genos.numberOfTaxa(); t++) {
			count_1 = 0;
			count_0 = 0;
			for (int site = 0; site < genos.numberOfSites(); site++) {
				int ct = crossType.get(site);
				byte a0 = 0; // less significant //minor allele
				byte a1 = 0; // most significant //major allele
				String g = null;
				String h = null;
				String gene = genos.genotypeAsString(t, site);
				if (ct == 1 || ct == 2 || ct == 3) {
					a1 = exp_f1[site][0];
					a0 = exp_f1[site][1];

					g = NucleotideAlignmentConstants.getNucleotideIUPAC(GenotypeTableUtils.getDiploidValue(a1, a0));

					if (g == gene)
						count_1++;
					else
						count_0++;
					continue;
				} else if (ct == 4 || ct == 5 || ct == 6) {
					a1 = exp_f1[site][0];
					a0 = exp_f1[site][1];

					g = NucleotideAlignmentConstants.getNucleotideIUPAC(a1);
					h = NucleotideAlignmentConstants.getNucleotideIUPAC(a0);

					if (gene == g || gene == h)
						count_1++;
					else
						count_0++;
					continue;
				}
			}
//			System.out.print((count_1 > count_0) ? "1 " : "0 ");
			if (count_1 > count_0) {
				scores_auto[t] = 1;
			} else
				scores_auto[t] = 0;
		}
		return scores_auto;
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
//		System.out.println("\nTook " + (endTime - startTime) + " ns");
	}

	private static void start(String hmpFile, String pedFile, double cutOff) throws InterruptedException {
//		System.out.println("Cut-off: " + cutOff);

		try {
			// READ FILE
//			Step1.listOfGenos = ReadFile(hmpFile, pedFile);
			// GET KEYS FROM CORRESPONDING PEDIGREE FILE
//			ArrayList<String> keys = PedigreeFileInfo.getKeys();

//			System.out.println(listOfGenos.size());

//			for(int i = 0; i < listOfGenos.size(); i++) {
//				genos = listOfGenos.get(i);
//				System.out.println(genos.hasGenotype());
//				ArrayList<Integer> crossType = Step1.GetCrossType(genos);
//				byte[][][] solution = Step1.InferGenotype(genos, crossType);
//				//TEST SOLUTIONS: f1ExpPerCross for subroutine 3.2a
//				testSolution(genos, solution[0], crossType);
//				ProgramEnded();
//				break;
//}

			// IMPLEMENT STEP 1
			System.out.print("In " + hmpFile + ": ");
			genos = Step1.ReadHMPFile(hmpFile);
//			ArrayList<String> keys = PedigreeFileInfo.getKeys();

//			for(int i = 0; i < keys.size(); i++) {
//				System.out.println(keys.get(i));
//			}
			ArrayList<Integer> crossType = Step1.GetCrossType(genos);
			if(crossType == null) { //no good markers were present
//				System.out.println("(MARKER_ERR) No good markers present. Handling will be done soon. \n");
//				Step1.scores_auto = {}
				Step1.scores_auto = new int[genos.numberOfTaxa()];
				for(int i = 0; i < genos.numberOfTaxa(); i++)
					Step1.scores_auto[i] = 0;
//				System.exit(0);
			}else {
			byte[][][] solution = Step1.InferGenotype(genos, crossType);
//			byte[][] f1ExpPerCross = solution[0]; //FOR SUBROUTINE 3.2a
//			byte[][] selfExpPerCross = solution[1]; //for subroutine 3.2b	

			// TEST SOLUTIONS: f1ExpPerCross for subroutine 3.2a
			Step1.scores_auto = testSolution(genos, solution[0], crossType); //f1's
			//SELFING LINES
//			Step1.scores_auto_selfing = testSolution(genos, solution[1], crossType, 0); //selfing lines

//			System.out.println("\n");
//			for(int i = 0; i < crossType.size(); i++) System.out.println(crossType.get(i));						
			}
			ProgramEnded();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public int[] getAutomatedScores() {
		return Step1.scores_auto;
	}
	
	public int[] getSelfingScores() {
		return Step1.scores_auto_selfing;
	}
}
