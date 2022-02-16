int getFlavorBTV(int flavor){
    //Maps hadronFlavor to BTV flavor:
    //Note the flavor convention: hadronFlavor is b = 5, c = 4, f = 0
    //Convert them to the btagging group convention of 0, 1, 2
    int flavor_btv;
    if (abs(flavor) == 5) flavor_btv = 0;
    else if (abs(flavor) == 4) flavor_btv = 1;
    else if (abs(flavor) == 0 || abs(flavor) == 1 || abs(flavor) == 2 || abs(flavor) == 3 || abs(flavor) == 21) flavor_btv = 2;
    else{
        cout<<"WARNING: Unknown flavor "<<flavor<<"setting b-tagging SF to -1!"<<endl;
        return -1;
        }
    return flavor_btv;
}

RVec<RVec<float>> btagSF(string sample, string WP, rvec_f Jet_pt, rvec_f Jet_eta, rvec_i Jet_hadronFlavour, rvec_f Jet_btagDeepFlavB, string era, string wp){
    RVec<RVec<float>> result;
    float max_abs_eta = 2.4;
    string discr;
    //if (algo == "csvv2") discr = "btagCSVV2";
    //else if (algo == "deepcsv") discr = "btagDeepB";
    //else if (algo == "cmva") discr = "btagCMVA";
    //else if (algo == "deepjet") discr = "btagDeepFlavB";
    //else cout<<"ERROR: Invalid algorithm "<< algo<<endl;
    
    BTagCalibrationReader reader;
    if (era == "Legacy2016"){
        if (wp == "L") reader = reader_0_Legacy2016;
        else if (wp == "M") reader = reader_1_Legacy2016;
        else reader = reader_2_Legacy2016;
    }
    else if (era == "2017"){
        if (wp == "L") reader = reader_0_2017;
        else if (wp == "M") reader = reader_1_2017;
        else reader = reader_2_2017;
    }
    else if (era == "2018"){
        if (wp == "L") reader = reader_0_2018;
        else if (wp == "M") reader = reader_1_2018;
        else reader = reader_2_2018;
    }
    else if (era == "UL2016APV"){
        if (wp == "L") reader = reader_0_UL2016APV;
        else if (wp == "M") reader = reader_1_UL2016APV;
        else reader = reader_2_UL2016APV;
    }
    else if (era == "UL2016"){
        if (wp == "L") reader = reader_0_UL2016;
        else if (wp == "M") reader = reader_1_UL2016;
        else reader = reader_2_UL2016;
    }
    else if (era == "UL2017"){
        if (wp == "L") reader = reader_0_UL2017;
        else if (wp == "M") reader = reader_1_UL2017;
        else reader = reader_2_UL2017;
    }
    else {
        if (wp == "L") reader = reader_0_UL2018;
        else if (wp == "M") reader = reader_1_UL2018;
        else reader = reader_2_UL2018;
    }
    
    for(int i = 0; i<Jet_pt.size(); i++){
        RVec<float> sfs;
        float pt = Jet_pt[i];
        float eta = Jet_eta[i];
        int flavor_btv_int = getFlavorBTV(Jet_hadronFlavour[i]);
        
        float discr = Jet_btagDeepFlavB[i];
        
        float epsilon = 1.e-3;
        
        if (eta <= -max_abs_eta) eta = -max_abs_eta + epsilon;
        if (eta >= +max_abs_eta) eta = +max_abs_eta - epsilon;

        if(flavor_btv == 0){
            float sf = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt);
            float sf_up = reader.eval_auto_bounds("up", BTagEntry::FLAV_B, eta, pt);
            float sf_down = reader.eval_auto_bounds("down", BTagEntry::FLAV_B, eta, pt);
        }
        else if(flavor_btv == 1){
            float sf = reader.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt);
            float sf_up = reader.eval_auto_bounds("up", BTagEntry::FLAV_C, eta, pt);
            float sf_down = reader.eval_auto_bounds("down", BTagEntry::FLAV_C, eta, pt);
        }
        else {
            float sf = reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, pt);
            float sf_up = reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG, eta, pt);
            float sf_down = reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG, eta, pt);
        }
        
        float sf, sf_up, sf_down;
        
        // check if SF is OK
        if (sf < 0.01) sf = 1.;
        if (sf_up < 0.01) sf_up = 1.;
        if (sf_down < 0.01) sf_down = 1.;
        
        sfs[0] = sf;
        sfs[1] = sf_up;
        sfs[2] = sf_down;
        
        result.emplace_back(sfs);
    }
    
    return result;
}
