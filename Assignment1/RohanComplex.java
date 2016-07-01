import edu.gwu.lintool.ComplexNumber;

/**
 * Rohan Sunder
 * CSCI 6342 - Linear Algebra: A Computational Approach
 * Assignment 1
 * @author rsunder
 *
 */
public class RohanComplex extends ComplexNumber {

	
	public RohanComplex(double a, double b) {
		this.re = a;
		this.im = b;
	}
	
	//Empty constructor
	public RohanComplex() {
		this.re = 0.0;
		this.im = 0.0;
	}
	
	@Override
	public ComplexNumber add(ComplexNumber arg0) {
		
		RohanComplex c = (RohanComplex) arg0;
		
		c.re = this.re + arg0.re;
		c.im = this.im + arg0.im;
		
		return c;
	}

	@Override
	public double angle() {
		
		if (this.re == 0) {
			return 0;
		}
		
		this.angle = Math.atan2(this.im, this.re);
		if (this.angle < 0) {
			this.angle = this.angle + 2 * (Math.PI);
		}
		
		return (this.angle);
	}

	@Override
	public ComplexNumber conjugate() {
		
		RohanComplex c = new RohanComplex();
		
		c.re = this.re;
		c.im = this.im * -1;
		
		return c;
	}

	@Override
	public double magnitude() {
		
		return Math.sqrt((this.re * this.re) + (this.im * this.im));
		
	}

	@Override
	public ComplexNumber mult(ComplexNumber arg0) {
		
		double first, inner, outer, last;
		RohanComplex c  = (RohanComplex)arg0;
		
		//Use Foil
		first = this.re * arg0.re;
		outer = this.re * arg0.im;
		inner = this.im * arg0.re;
		last = this.im * arg0.im * -1;
		
		c.re = first + last;
		c.im = outer + inner;
		
		return c;
	}

	@Override
	public ComplexNumber mult(double arg0) {
		
		RohanComplex c = new RohanComplex();
		
		c.angle = this.angle;
		
		c.re = arg0 * this.re;
		c.im = arg0 * this.im;
		
		return c;
	}

	@Override
	public ComplexNumber pow(int arg0) {
		
		RohanComplex c = new RohanComplex();
		System.out.println(arg0 + "");
		
		c.re = this.re;
		c.im = this.im;

		for (int i = 0; i < arg0-1; i++) {
			c = (RohanComplex) this.mult(c);
		}
			
		return c;
		
	}

	@Override
	public ComplexNumber sub(ComplexNumber arg0) {
		
		RohanComplex c = (RohanComplex) arg0;
		
		c.re = this.re - arg0.re;
		c.im = this.im - arg0.im;
		
		return c;
	}

}
