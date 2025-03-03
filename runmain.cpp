#include "runmain.h"

using namespace std;


void readinput(int *args, int proc_id, int proc_size)
{
	ifstream file;
	string c;
	file.open("1_Input.txt");
	if(file==NULL)
	{
		cout<<"No Input file present!"<<endl;
		exit(0);
	}
	for(int j=0;j<10;j++)
	{
		getline(file,c,';');
		file>>args[j];
	}
	file.close();
	args[10]=proc_id;
	args[11]=proc_size;
}

void readplasmainput(double* args)
{

	ifstream file;
	string c;
	file.open("1_Plasma.txt");
	if(file==NULL)
	{
		cout<<"No Plasma input file present!"<<endl;
		exit(0);
	}
	for(int j=0;j<16;j++)
	{
		getline(file,c,';');
		file>>args[j];
	}
	file.close();
}

int main (int argc, char *argv[])
{
	int proc_id, proc_size;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_size);
	int coordinatedoutput=0;
	int inputs[12];
	double plasmain[16];

	for(int i=0;i<proc_size;i++)
	{
		if(proc_id==i)
		{
			readinput(inputs,proc_id,proc_size);
			readplasmainput(plasmain);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	ofstream file,Fp1,Fp2;

	if(proc_id==0)
	{

		cout<<"New Program has started"<<endl;
		//---------------------------------------------------------------------------------------------------------------------------------------------------------------
		char q='"';

		if(inputs[7]==1)
		{
			file.open("2_Time_characteristics.dat");

			file<<"Variables="<<q
					<<"T"<<q<<","<<q
					<<"<greek>D</greek>V"<<q<<","<<q
					<<"j<sub>ne<sub>drift</sub></sub>"<<q<<","<<q
					<<"j<sub>ne<sub>drift+diff</sub></sub>"<<q<<","<<q
					<<"j<sub>ne<sub>drift+diff</sub></sub>+j<sub>displacement</sub>"<<q<<","<<q
					<<"R<sub>ni</sub>"<<q<<","<<q
					<<"R<sub>ne</sub>"<<q<<","<<q
					<<"R<sub>phi</sub>"<<q<<","<<q
					<<"<greek>D</greek>t<sub>in_use</sub>"<<q<<","<<q
					<<"<greek>D</greek>t<sub>diffusion</sub>"<<q<<","<<q
					<<"<greek>D</greek>t<sub>pos_presv</sub>"<<q<<","<<q
					<<"<greek>D</greek>t<sub>plasma</sub>"<<q<<","<<q
					<<endl;

			Fp1.open("2_Potdrop_Non-dimensional.dat");
			Fp1<<"Variables="<<q
					<<"X"<<q<<","<<q
					<<"e<greek>D</greek><greek>f</greek>/(k<sub>B</sub>T<sub>e</sub>)"<<q<<","<<q<<endl;

			Fp2.open("2_Potdrop_potential.dat");
			Fp2<<"Variables="<<q
					<<"X"<<q<<","<<q
					<<"<greek>f</greek>"<<q<<","<<q<<endl;

		}
		else
		{
			file.open("2_Time_characteristics.dat",fstream::app);
			Fp1.open("2_Potdrop_Non-dimensional.dat",fstream::app);
			Fp2.open("2_Potdrop_potential.dat",fstream::app);
		}
		//---------------------------------------------------------------------------------------------------------------------------------------------------------------
	}



	globalvar gv;
	gv.allocate(inputs);


	double www[10];
	for(int j=0;j<10;j++)www[j]=0.0;
	int outputfreq=pow(10,7);
	//int outputfreq=2*pow(10,5);
	int raiseflag=0;
	double conv_crit=plasmain[15];
	double t,tmp;
	double T,DT;
	double R[6];
	int count=0;

	gv.initialization(&Fp1,&Fp2,plasmain);
	gv.matrix_assemble();
	MPI_Barrier(MPI_COMM_WORLD);
	gv.DTreset();
	DT=gv.DT();
	T=gv.T();
	t=gv.starttime();
	tmp=gv.TMP();

	while(t<=T)
	{
		www[0]=clock();

		if(count%500000==0 && proc_id==0)
		{
			cout<<"Iteration : "<<count<<" , t = "<<t;
			cout<<", Timestep in use = "<<DT;
			cout<<", Potential drop = "<<gv.retV();
			cout<<", solve time : "<<www[1]/CLOCKS_PER_SEC<<"s per iteration"<<endl;
			//cout<<"Time steps: convection = "<<gv.DTconvection()<<", diffusion: "<<gv.DTdiffusion()<<", pospresv: "<<gv.DTpospresv()<<", plasma: "<<gv.DTplasma()<<endl;
		}

		count++;
		gv.poisson_solver();
		gv.solve();
		//---------------------------------------------------------------------------------------------------------------------------------------------------------------
		
		if(count%1000==0 && count>10)
		{
			MPI_Status stat;
			gv.residuals(R);
			if(proc_id>0)
				MPI_Send(&R[0],6,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
			if(proc_id==0)
			{
				double s[6],sum[6];for(int j=0;j<6;j++) sum[j]=0.0;
				for(int op=1;op<inputs[9];op++)
				{
					MPI_Recv(&s[0],6,MPI_DOUBLE,op,1,MPI_COMM_WORLD,&stat);
					for(int j=0;j<6;j++) sum[j]=sum[j]+s[j];
				}
				for(int j=0;j<6;j++) sum[j]=(sum[j]+R[j]);
				for(int j=0;j<3;j++) sum[j]=pow((2.0*inputs[0]-2.0)*sum[j],0.5)/sum[j+3];
				if(sum[0]<conv_crit && sum[1]<conv_crit && sum[2]<conv_crit)
				{
					if(!(count<10000 && inputs[7]==2))
					{
						cout<<"Solution has converged at time: "<<t+gv.DT()<<endl;
						MPI_Abort(MPI_COMM_WORLD,1);
					}
				}

				file<<t<<" "<<gv.retV()<<" "<<gv.ret_current()<<" "<<gv.ret_conduction_current()<<" "<<gv.ret_total_current()<<" "
						<<sum[0]<<" "<<sum[1]<<" "<<sum[2]<<" "
						<<gv.DT()<<" "<<gv.DTdiffusion()<<" "<<gv.DTpospresv()<<" "<<gv.DTplasma()
						<<endl;
			}
		}
		
		//---------------------------------------------------------------------------------------------------------------------------------------------------------------

		gv.update();
		gv.updatet();

		DT=gv.DT();
		t=t+DT;
		count++;
		gv.DTreset();


		if(count%outputfreq/10==0)
			gv.back_up();


		if(count%outputfreq==0 || count==pow(10,coordinatedoutput))
		{
			gv.output();
			coordinatedoutput++;
			if(t>5*pow(10,-5) && raiseflag==0)
			{
				outputfreq=outputfreq*10;
				raiseflag=1;
			}
		}

		www[1]=((count-1)*www[1]+clock()-www[0])/(count);


		MPI_Barrier(MPI_COMM_WORLD);

	}//end of time stepping
	gv.output();
	gv.back_up();

	gv.deallocate();

	file.close();
	Fp1.close();
	Fp2.close();

	cout<<"Program has ended"<<endl;

	MPI_Finalize();

	return 0;
}








