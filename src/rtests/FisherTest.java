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

public class FisherTest {
    public static void main(String[] args) {
	if (!Rengine.versionCheck()) {
	    System.err.println("** Version mismatch - Java files don't match library version.");
	    System.exit(1);
	}
	Rengine re=new Rengine(args, false, null);
		// the engine creates R is a new thread, so we should wait until it's ready
        if (!re.waitForR()) {
            System.out.println("Cannot load R");
            return;
        }
		try {
			String sample_m = "matrix(c(3, 1, 1, 3), nrow = 2, dimnames = list(Guess = c(\"Milk\", \"Tea\"),Truth = c(\"Milk\", \"Tea\")))";
			re.eval("TeaTasting=" + sample_m);
			re.eval("result=fisher.test(TeaTasting, alternative = \"greater\")");
			Double pval = re.eval("result$p.value").asDouble();
			System.out.println("P-value is = " + pval);
		}catch(Exception e) {}
    }
}
