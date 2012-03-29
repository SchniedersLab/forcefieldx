package ffx.algorithms;

import ffx.numerics.Potential;

public class Respa extends Integrator {
    private double x[];
    private double v[];
    private double a[];
    private double a_alt[];
    private double aPrevious[];
    private double mass[];
    private int nVariables;
    private double dt;
    private double dt2_8;
    private double dt_8;
    private double dt_2;
    private double dalt;
    private double dta;
    private double dta_2;
    private int nalt;
    private double eps = .00000001;
    public double halfStepEnergy = 0;
    
    public Respa(int nVariables, double x[], double v[], double a[],
                        double aPrevious[], double mass[]) {
        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.a = a;
        this.aPrevious = aPrevious;
        this.mass = mass;
        a_alt = new double[nVariables];
        dt = 1.0;
        dt_2 = .5 * dt;
        dalt = .00025;
        nalt = ((int) (dt/(dalt+eps))) + 1;
        dalt = (double) nalt;
        dta = dt/dalt;
        dta_2 = .5 * dta;
    }

	@Override
	public void halfStep(Potential potential) {
        double gradient[] = new double[nVariables];
		for (int i = 0; i < nVariables; i++) {
            v[i] += a[i] * dt_2;
        }

        for(int j = 0; j < nalt; j++){
    		for (int i = 0; i < nVariables; i++) {
    			x[i] += v[i] * dta_2;
            }
    		potential.turnFastOnSlowOff(1);
        	halfStepEnergy = potential.energyAndGradient(x, gradient);
            for (int i = 0; i < nVariables; i++) {
                a_alt[i] = -Thermostat.convert * gradient[i] / mass[i];
                v[i] += a_alt[i] * dta;
                x[i] += v[i] * dta_2;
            }
        }
        potential.turnFastOnSlowOff(0);
	}
	@Override
	public void fullStep(double[] gradient) {
        for (int i = 0; i < nVariables; i++) {
            a[i] = -Thermostat.convert * gradient[i] / mass[i];
            v[i] += a[i] * dt_2;
        }
	}
	
	@Override
	public void setTimeStep(double dt) {
        this.dt = dt;
        dt_2 = .5 * dt;
        dalt = .00025;
        nalt = ((int) (dt/(dalt+eps))) + 1;
        dalt = (double) nalt;
        dta = dt/dalt;
        dta_2 = .5 * dta;
	}

}

