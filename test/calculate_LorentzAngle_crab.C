{
	setTDRStyle();
	
	TFile *f = new TFile("lorentzangle_Drell_Yan_10pb_490000ev.root", "READ");
	f->cd();
	
	TF1 *f1 = new TF1("f1","[0] + [1]*x",50., 235.); 
	f1->SetParName(0,"p0");
	f1->SetParName(1,"p1");
	f1->SetParameter(0,0);
	f1->SetParameter(1,0.4);
	ofstream fLorentzFit( "lorentzFit.txt", ios::trunc );
	fLorentzFit.precision( 4 );
	fLorentzFit << "module" << "\t" << "layer" << "\t" << "offset" << "\t" << "error" << "\t" << "slope" << "\t" << "error" << "\t" "rel.err" << "\t" "pull" << "\t" << "chi2" << "\t" << "prob" << endl;
	//loop over modlues and layers to fit the lorentz angle
	for( int i_layer = 1; i_layer<=3; i_layer++){
		for(int i_module = 1; i_module<=8; i_module++){
			TString name = "h_drift_depth_adc_layer";
			name += i_layer;
			name += "_module";
			name += i_module;
			TH2F * h_drift_depth_adc = (TH2F*) f->Get(name);
			TString name = "h_drift_depth_adc2_layer";
			name += i_layer;
			name += "_module";
			name += i_module; 
			TH2F * h_drift_depth_adc2 = (TH2F*) f->Get(name);
			TString name = "h_drift_depth_noadc_layer";
			name += i_layer;
			name += "_module";
			name += i_module;
			TH2F * h_drift_depth_noadc = (TH2F*) f->Get(name);
			int hist_drift_ = h_drift_depth_adc->GetNbinsX();
			int hist_depth_ = h_drift_depth_adc->GetNbinsY();
			float min_depth_ = h_drift_depth_adc->GetYaxis()->GetXmin();
			float max_depth_ = h_drift_depth_adc->GetYaxis()->GetXmax();
			float min_drift_ = h_drift_depth_adc->GetXaxis()->GetXmin();
			float max_drift_ = h_drift_depth_adc->GetXaxis()->GetXmax();
			TH1F * h_mean = new TH1F("h_mean","h_mean", hist_depth_, min_depth_, max_depth_);
			TH1F * h_drift_depth_adc_slice_ = new TH1F("h_slice","h_slice", hist_drift_, min_drift_, max_drift_);
			//loop over bins in depth (z-local-coordinate) (in order to fit slices)
			for( int i = 1; i <= hist_depth_; i++){
// 				findMean(i, (i_module + (i_layer - 1) * 8));
				double nentries = 0;
	
				h_drift_depth_adc_slice_->Reset("ICE");
				
				// determine sigma and sigma^2 of the adc counts and average adc counts
					//loop over bins in drift width
				for( int j = 1; j<= hist_drift_; j++){
					if(h_drift_depth_noadc->GetBinContent(j, i) >= 1){
						double adc_error2 = (h_drift_depth_adc2->GetBinContent(j,i) - h_drift_depth_adc->GetBinContent(j,i)*h_drift_depth_adc->GetBinContent(j, i) / h_drift_depth_noadc->GetBinContent(j, i)) /  h_drift_depth_noadc->GetBinContent(j, i);
						h_drift_depth_adc_slice_->SetBinContent(j, h_drift_depth_adc->GetBinContent(j,i));
						h_drift_depth_adc_slice_->SetBinError(j, sqrt(adc_error2));
						nentries += h_drift_depth_noadc->GetBinContent(j,i);	
					}else{
					 	h_drift_depth_adc_slice_->SetBinContent(j, h_drift_depth_adc->GetBinContent(j,i));
						h_drift_depth_adc_slice_->SetBinError(j, 0);
					}
				} // end loop over bins in drift width
					
				double mean = h_drift_depth_adc_slice_->GetMean(1); 
				double error = 0;
				if(nentries != 0){
					error = h_drift_depth_adc_slice_->GetRMS(1) / sqrt(nentries);
				}
					
				h_mean->SetBinContent(i, mean);
				h_mean->SetBinError(i, error);	
			}// end loop over bins in depth 
			h_mean->Fit(f1,"ERQ");
			double p0 = f1->GetParameter(0);
			double e0 = f1->GetParError(0);
			double p1 = f1->GetParameter(1);
			double e1 = f1->GetParError(1);
			double chi2 = f1->GetChisquare();
			double prob = f1->GetProb();
			delete h_mean;
			delete h_drift_depth_adc_slice_;
			fLorentzFit << i_module << "\t" << i_layer << "\t" << p0 << "\t" << e0 << "\t" << p1 << "\t" << e1 << "\t" << e1 / p1 *100. << "\t"<< (p1 - 0.424) / e1 << "\t"<< chi2 << "\t" << prob << endl;
		}
	}
	
}
