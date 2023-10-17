void GetData()
{
TFile *f = new TFile("EffArea_Signal_0_0.root");

TH1D *hback = (TH1D*)f->Get("EffectiveArea_5");


ofstream file_out;
file_out.open("EffArea.dat");
//file_out=open("DiffuseGamma.txt");

for (Int_t i=1; i<hback->GetEntries(); i++) {
    cout<<hback->GetBinCenter(i)<<" "<<hback->GetBinLowEdge(i)<<" "<<hback->GetBinLowEdge(i+1)<<" "<<hback->GetBinContent(i)<<endl;
     file_out <<hback->GetBinCenter(i)<<" "<<hback->GetBinLowEdge(i)<<" "<<hback->GetBinLowEdge(i+1)<<" "<<hback->GetBinContent(i)<<endl;

  }

}
