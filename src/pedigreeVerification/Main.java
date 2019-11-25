package pedigreeVerification;

import java.util.ArrayList;

public class Main{
	public static void main(String[] args) {
		PreprocessHMP n = new PreprocessHMP(args);
		if(n.getVal()) {
			Step1 method1 = new Step1(args);
			int[] ascore = method1.getAutomatedScores();
			ArrayList<Integer> mscore = n.getManualScores();
			//IDENTIFY TRUE AND FALSE NEGATIVES!!!!!
//			int tp=0, tn=0, fp=0, fn=0;
			
			float count_1 = 0;
			float count_0 = 0;
			System.out.print("\nmanual: ");
			for(int i = 0; i < mscore.size(); i++) System.out.print(mscore.get(i)+" ");
			System.out.print("\nauto:   ");
			for(int i = 0; i < mscore.size(); i++) System.out.print(ascore[i]+" ");
			int tp=0;
			int tn=0;
			int fp=0; 
			int fn=0;
			System.out.println("\n\nConsidered True F1's (upon automation): ");
			for(int i = 0; i < ascore.length; i++) {
				if(mscore.get(i) == 1 && ascore[i] == 1) {
					tp++;
					System.out.println((n.getFilename().split(".hmp.txt")[0])+"-"+(i+1));
				}
				else if(mscore.get(i) == 0 && ascore[i] == 0)	tn++;
				else if(mscore.get(i) == 0 && ascore[i] == 1)	fp++;
				else if(mscore.get(i) == 1 && ascore[i] == 0)	fn++;
			}
			System.out.print("\nTP:"+tp+" TN:"+ tn+" FP:"+fp+" FN:"+ fn+"\n\n");
			System.exit(0);
//			//count yung mga number 1
//			int count_ng_1 = 0;
//			int count_ng_0 = 0;
//			for(int j = 0; j < mscore.size(); j++) {
//				if(mscore.get(j) == 1) count_ng_1++; 
//				else count_ng_0++;
//			}
//			Double rate = (double) (count_1 / count_ng_1);
//			Double rate0 = (double) (count_0 / count_ng_0);
//			System.out.println((int) (rate * 100 ) + "% (F1 ONLY)\n"); //IMPORTANT: INDICATE IN LOG FILE THAT ONLY F1'S ARE CONSIDERED. SELFING LINES ARE NOT COVERED BY THE SCOPE OF THE SPECIAL PROBLEM.
		}else {
			if(n.getFilename() != null && !(n.getFilename().equals("whole_snp_seq.hmp.txt"))) 
				System.out.println("!! For file " + n.getFilename() + ": (DATA_ERR) Data anomalies were encountered. No processing was done.\n");
		}
	}
}