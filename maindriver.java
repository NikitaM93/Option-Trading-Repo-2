import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.lang.Integer;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import java.awt.GridLayout;
import javax.swing.JLabel;
import java.awt.FlowLayout;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JTextField;
import javax.swing.JPasswordField;
import javax.swing.JOptionPane;

public class maindriver{
	public ArrayList<MarketDate> dates=new ArrayList<MarketDate>();
	public ArrayList<Double> logreturns=new ArrayList<Double>();
	public ArrayList<Double> averagelogs=new ArrayList<Double>();
	public ArrayList<Double> variancecalcs=new ArrayList<Double>();
	public double annualvolprime;
	public double pi=3.14159265358979;
	public double annualvolprimeimplied;
	public long daysToMaturity=timeToMaturity(2014,4,25,0,1);
	public int daysToMaturityPrime2=(int) daysToMaturity;
	
    
	
	public static void main(String[] args){
		maindriver md = new maindriver(/* constructor args here */);
		String temp="";

		md.reader();
		
		
		//MAKES VOL, vol is saved in annualvolprime variable. 
    	md.logreturn(dates);
    	md.averagelog(logreturns);
    	md.variancecalc(logreturns);
    	md.varianceaverage(variancecalcs);
    	md.stdev(varianceaverage(variancecalcs));
    	md.annualvol(stdev(varianceaverage(variancecalcs)));
  





   // 	System.out.println("annualvolprime>>>>>>>>>>>>>>>>>>>>>>>>>>>"+annualvolprime);
        
    			//************************************************************************************************************************************
    			//***********Test of Wiki BS**********************************************************************************************************
    //    		long daysToMaturityPrime=timeToMaturity(2014,4,25,0,1);
    //    		long call=europeanCallPrice(annualvolprime,daysToMaturityPrime,539.19,525.00,0.001);
    			//System.out.println(annualvolprime);
    			//System.out.println(daysToMaturity);
    //			long put=europeanPutPrice(annualvolprime,daysToMaturityPrime,539.19,525.00,0.001);
    //			System.out.println("Price of E.Call: $"+call);
    //			System.out.println("Price of E.Put: $"+put);
    			//************************************************************************************************************************************
    	
    	System.out.println("--------------------------------\n---------------SG's No Beta-----\n---------------Hedge Alpha------\n--------------------------------\n");
    	//OptionPrice(spot,strike,daystoMaturity,vol,rate,q,optiontype as in c or p)
    	//long daysToMaturity=timeToMaturity(2014,4,25,0,1);
    	String c="c";
    	String p="p";
    	//int daysToMaturityPrime2=(int) daysToMaturity;
    	
    	double SgOption=md.OptionPrice(531.17,510.00,daysToMaturityPrime2,annualvolprime,0.0012,0,c);
    	//long SgOptionPut=OptionPrice(64.82,59.50,daysToMaturityPrime2,annualvolprime,0.0012,0,p);
    	
    	double SgGreek=md.greeks(544.70,525.00,daysToMaturityPrime2,annualvolprime,0.0013,0.0001,c);
    	System.out.println("SG's price of Call: $"+SgOption);
    	//System.out.println("SG's price of Put: $"+SgOptionPut);
    	//System.out.println("TIMMMMMMMMMMMME"+daysToMaturityPrime2);
    	
    	
    	double impliedvol=md.annualvolimplied(SgOption, 24.10);
    	System.out.println("SG's implied annualvol for Call: "+impliedvol);
    	
    	
    	
    }
	public double OptionPrice(double spot, double strike, double NbExp, double vol, double rate, double q, String optionType){
		double v, d1, d2, Nd1, Nd2, T, result;
		//System.out.println(spot);
		//System.out.println(strike);
		//System.out.println(NbExp);
		//System.out.println(vol);
		//System.out.println(rate);
		//System.out.println(q);
		//System.out.println(optionType);
		
		if(NbExp<0)
			return 0;
		
		T=NbExp/365;
		if(NbExp==0){
			if(optionType.equals("c")){
				System.out.println("TESTTESTTESTTEST"+(long) (Math.max(spot-strike,0)));
				return (double) (Math.max(spot-strike,0));
			}
			else{
				System.out.println("TESTTESTTESTTEST"+(long) (Math.max(spot-strike,0)));
				return (double) (Math.max(strike-spot,0));	
			}
		}
		
		d1=((Math.log(spot/strike))+(rate-q+(vol*vol)/2)*T)/(vol*Math.sqrt(T));
		//System.out.println("TESTTESTTESTTEST");
		//System.out.println((Math.log(spot/strike)));
		//System.out.println(vol*Math.sqrt(T));
		//System.out.println("TESTTESTTESTTEST"+d1);
		
		d2=d1-vol*Math.sqrt(T);
		//System.out.println("TESTTESTTESTTEST2  "+d2);
		Nd1=cdnf(d1);
		//System.out.println("TESTTESTTESTTEST Nd1  "+Nd1);
		Nd2=cdnf(d2);
		//System.out.println("TESTTESTTESTTEST Nd2  "+Nd2);
		if(optionType.equals("c")){ //call option
			//System.out.println((long) (spot*Math.exp(-q*T)*Nd1-strike*Math.exp(-rate*T)*Nd2));
			return (double) (spot*Math.exp(-q*T)*Nd1-strike*Math.exp(-rate*T)*Nd2);
		}
		
		else //put option
			return (double) (-spot*Math.exp(-q*T)*(1-Nd1)+strike*Math.exp(-rate*T)*(1-Nd2));
	}
	
	public double greeks(double spot, double strike, double NbExp, double vol, double rate, double q, String optionType){
		double v, d1, d2, Nd1, Nd2, T, result, dS, dv, dr;
		double dt;
		double delta, gamma, vega, theta, rho;
		
		dS=0.01; //0.01 point move in spot
		dv=0.0001; //0.01% move in vol
		dt=1; //1 day
		dr=0.0001; //1bps move
		
		if(NbExp<0)
			//System.out.println("TESTTESTTESTTEST");
			return 0;
		
		double x=(double)(OptionPrice(spot+dS,strike,NbExp,vol,rate,q,optionType));
		double x2=(double)OptionPrice(spot-dS,strike,NbExp,vol,rate,q,optionType);
		//System.out.println(x);
		//System.out.println(x2);
		//System.out.println(dS);
		//System.out.println((x-x2)/(2*dS));
		
		delta=(double) ((OptionPrice(spot+dS,strike,NbExp,vol,rate,q,optionType)-OptionPrice(spot-dS,strike,NbExp,vol,rate,q,optionType))/(2*dS));
		
		gamma=(double) ((OptionPrice(spot+dS,strike,NbExp,vol,rate,q,optionType)-2*OptionPrice(spot,strike,NbExp,vol,rate,q,optionType)+OptionPrice(spot-dS,strike,NbExp,vol,rate,q,optionType))/(dS*dS));
		vega=(double) ((OptionPrice(spot,strike,NbExp,vol+dv,rate,q,optionType)-OptionPrice(spot,strike,NbExp,vol-dv,rate,q,optionType))/(2*dv)/100);
		rho=(double) ((OptionPrice(spot,strike,NbExp,vol,rate+dr,q,optionType)-OptionPrice(spot,strike,NbExp,vol,rate-dr,q,optionType))/(2*dr)/1000);

		if(NbExp==0){
			//System.out.println("TESTTESTTESTTEST");
			theta=0;
		}
		else 
			theta=(double) ((OptionPrice(spot,strike,NbExp-dt,vol,rate,q,optionType)-OptionPrice(spot,strike,NbExp+dt,vol,rate,q,optionType))/(2*dt));
		
		System.out.println("delta>>>>>>"+delta);
		System.out.println("gamma>>>>>>"+gamma);
		System.out.println("vega>>>>>>"+vega);
		System.out.println("rho>>>>>>"+rho);
		return 0;
	}
	
	public long timeToMaturity(int year, int month, int day, int hour, int minute){
		Date maturity=new GregorianCalendar(year, month, day, hour, minute).getTime();
    	/** Today's date */
		Date today=new Date();
		// Get msec from each, and subtract.
		long diff=maturity.getTime()-today.getTime();
		long timeToMaturity=diff/ (1000 * 60 * 60 * 24);
		return timeToMaturity;
	}
	
	//Reference:http://stackoverflow.com/questions/442758/which-java-library-computes-the-cumulative-standard-normal-distribution-function
	public double cdnf(double x){
	    int neg=(x<0d)?1:0;
	    if (neg==1) 
	        x*=-1d;
	    double k=(1d/(1d+0.2316419*x));
	    double y=((((1.330274429*k-1.821255978)*k+1.781477937)*k-0.356563782)*k+0.319381530)*k;
	    y=1.0-0.398942280401*Math.exp(-0.5*x*x)*y;
	    //System.out.println((1d-neg)*y+neg*(1d-y));
	    return (1d-neg)*y+neg*(1d-y);
	}
	
	//For Call Options
	public double annualvolimplied(double modelOption, double realOption){
		double volimplied=stdev(varianceaverage(variancecalcs));
		double annualvolimplied=0;
		double real=realOption;
		double model=modelOption;
		
		if(real==model){
			annualvolimplied=annualvol2(volimplied);
		}
		
		else if(real>model){
			do{
				volimplied=volimplied+0.00001;
				model=OptionPrice(531.17,510.00,daysToMaturityPrime2,annualvol2(volimplied),0.0012,0,"c");
				annualvolimplied=annualvol2(volimplied);
			}while(real>model);
		}
		else
			do{
				volimplied=volimplied-0.00001;
				model=OptionPrice(531.17,510.00,daysToMaturityPrime2,annualvol2(volimplied),0.0012,0,"c");
				annualvolimplied=annualvol2(volimplied);
			}while(real<model);
		
		return annualvolimplied;	
	}
	
	public double annualvol2(double stdev){
		double x=Math.sqrt(252);
		double annualvol2=x*stdev;
		//System.out.println("AnnualVolImplied:"+annualvol);
		//annualvolprimeimplied=annualvol2;
		return annualvol2;	
	}
	
	
	//*********************Volatility Calc*******************************************************************************************
	public void logreturn(ArrayList dates){
		for(int i=0; i<dates.size()-1; i++){
			MarketDate day2=(MarketDate) dates.get(i);
			MarketDate day1=(MarketDate) dates.get(i+1);
			double x= Math.log(day1.getAdjClose()/day2.getAdjClose());
			//System.out.println("Day :"+ i+ ":" + x);
			logreturns.add(x);
		}
	}
	
	public double averagelog(ArrayList logreturns){
		double sum=0;
		double counter=0;
		for(int i=0; i<logreturns.size()-1; i++){
			double y=(Double) logreturns.get(i);
			sum+=y;
			//System.out.println(sum);
			counter++;	
		}
		double x=sum/counter;
		//System.out.println(x);
		return x;
	}
    
	public void variancecalc(ArrayList logreturns){
		for(int i=0; i<logreturns.size()-1; i++){
			double y=(Double) logreturns.get(i);
			double yprime=y-averagelog(logreturns);
			variancecalcs.add(Math.pow(yprime,2));
		}
	}
	
	public double varianceaverage(ArrayList variancecalc){
		double sum=0;
		double counter=0;
		for(int i=0; i<variancecalc.size()-1; i++){
			double y=(Double) variancecalc.get(i);
			sum+=y;
			//System.out.println(sum);
			counter++;	
		}
		double x=sum/counter;
		System.out.println("******************");
		System.out.println("VarAvrg:"+x);
		return x;
	}
	
	public double stdev(double varianceAvrg){
		double x=Math.sqrt(varianceAvrg);
		System.out.println("Stdev:"+x);
		return x;	
	}
	
	public double annualvol(double stdev){
		double x=Math.sqrt(252);
		double annualvol=x*stdev;
		System.out.println("AnnualVol:"+annualvol);
		annualvolprime=annualvol;
		return annualvol;	
	}
	//************************************************************************************************************************************
	public void reader(File f){        
        try{
            Scanner sc=new Scanner(f);

            int y=0;
            while(sc.hasNextLine()&&y<61){
                String line=sc.nextLine();
                String[] info=line.split(",");
                
                
                String date=info[0];
                double open=Double.parseDouble(info[1]);
                double high=Double.parseDouble(info[2]);
                double low=Double.parseDouble(info[3]);
                double close=Double.parseDouble(info[4]);
                double volume=Double.parseDouble(info[5]);
                double adjclose=Double.parseDouble(info[6]);
                
                MarketDate d=new MarketDate(date, open, high, low, close, volume, adjclose);
                dates.add(d);
                y++;
            }

            for(MarketDate d: dates){
                System.out.println(d.toString());
            }

        }
        catch (FileNotFoundException e) {         
        	e.printStackTrace();
        }
    }
	public void reader(){        
        try{
        	File file = new File("volcode.txt");
            Scanner sc=new Scanner(file);
            
            int y=0;
            while(sc.hasNextLine()&&y<61){
                String line=sc.nextLine();
                String[] info=line.split(",");
                
                
                String date=info[0];
                double open=Double.parseDouble(info[1]);
                double high=Double.parseDouble(info[2]);
                double low=Double.parseDouble(info[3]);
                double close=Double.parseDouble(info[4]);
                double volume=Double.parseDouble(info[5]);
                double adjclose=Double.parseDouble(info[6]);
                
                MarketDate d=new MarketDate(date, open, high, low, close, volume, adjclose);
                dates.add(d);
                y++;
            }

            for(MarketDate d: dates){
                System.out.println(d.toString());
            }

        }
        catch (FileNotFoundException e) {         
        	e.printStackTrace();
        }
    }
	
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//********************************Wiki implementation*********************************************************************************
	//**********************************of Black Scholes**********************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
		
	public long europeanCallPrice(double vol, long maturity, double spot, double strike, double rate){
		long Nd1S=(long) (cdnf(d1(vol, maturity, spot, strike, rate))*spot);
		//System.out.println(Nd1S);
		long Nd2KerTt= (long) (cdnf(d2(vol, maturity, spot, strike, rate))*strike*Math.exp(-rate*maturity));
		long callprice=Nd1S-Nd2KerTt;
		return callprice;
	}
	
	public long europeanPutPrice(double vol, long maturity, double spot, double strike, double rate){
		long KerTt= (long) (strike*Math.exp(-rate*maturity));
		long putprice=(long) (KerTt-spot+europeanCallPrice(vol, maturity, spot, strike, rate));
		return putprice;
	}
	
	public double d1(double vol, long maturity, double spot, double strike, double rate){
		double d1;
		double a1=1/(vol*Math.sqrt(maturity));
		double a2=Math.log(spot/strike);
		double a3=(rate+((vol*vol)/2))*(maturity);
		d1=(a2+a3)*a1;
		return d1;
	}
	
	public double d2(double vol, long maturity, double spot, double strike, double rate){
		double d2;
		double a1=1/(vol*Math.sqrt(maturity));
		double a2=Math.log(spot/strike);
		double a3=(rate-((vol*vol)/2))*(maturity);
		d2=(a2+a3)*a1;
		return d2;
	}
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
	//************************************************************************************************************************************
}

class MarketDate{

    private String date;
    private double open;
    private double high;
    private double low;
    private double close;
    private double volume;
    private double adjclose;
    private double weight;
    
    public MarketDate(String date, double open, double high, double low, double close, double volume, double adjclose){
    	this.setDate(date);
    	this.open=open;
    	this.high=high;
    	this.low=low;
    	this.close=close;
    	this.volume=volume;
    	this.adjclose=adjclose;     
    	weight=1;
    }
    public String getDate(){
    	return date;
    }
    public void setDate(String date){
    	this.date=date;
    }
    public double getOpen(){
    	return open;
    }
    public double getHigh(){
    	return high;
    }
    public double getLow(){
    	return low;
    }
    public double getClose(){
    	return close;
    }
    public double getVolume(){
    	return volume;
    }
    public double getAdjClose(){
    	return adjclose;
    }
    public String toString(){
    	return this.date + " " + this.open + " " + this.high + " " + this.low + " " + this.close + " " + this.volume + " " + this.adjclose;
    }
}
