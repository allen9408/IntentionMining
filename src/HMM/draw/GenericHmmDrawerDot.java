package HMM.draw;

import HMM.DataType.Hmm;


/**
 * An HMM to <i>dot</i> file converter.  See
 * <url>http://www.research.att.com/sw/tools/graphviz/</url>
 * for more information on the <i>dot</i> tool.
 * <p>
 * The command <tt>dot -Tps -o &lt;outputfile&gt; &lt;inputfile&gt;</tt>
 * should produce a Postscript file describing an HMM.
 */
public class GenericHmmDrawerDot extends HmmDrawerDot<Hmm<?>> {
	public void setMinAij (double val) {
		minimumAij = val;
	}
	
	public void setMinPI(double val) {
		minimumPi = val;
	}
	
	public void recoverDefault() {
		minimumAij = 0.01;
		minimumPi = 0.01;
	}
}
