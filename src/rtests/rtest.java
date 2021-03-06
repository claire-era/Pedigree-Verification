package rtests;
import java.io.*;
import java.awt.Frame;
import java.awt.FileDialog;

import java.util.Enumeration;

import org.rosuda.JRI.Rengine;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RList;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.RMainLoopCallbacks;

class rtest {
	public rtest() {}
	Rengine re=new Rengine(null, false, null);
	public double compute(){
		Double pval = null;
		try {
			String sample_m = "matrix(c(3, 1, 1, 3), nrow = 2, dimnames = list(Guess = c(\"Milk\", \"Tea\"),Truth = c(\"Milk\", \"Tea\")))";
			re.eval("TeaTasting=" + sample_m);
			re.eval("result=fisher.test(TeaTasting, alternative = \"greater\")");
			pval = re.eval("result$p.value").asDouble();
			System.out.println("P-value is = " + pval);
			System.out.println(pval);
			re.end();
			System.out.println(re.isAlive());
		} catch (Exception e) {
		}
       re.end();
       return pval;
    }
}
