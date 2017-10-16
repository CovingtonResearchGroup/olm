import sys,os
from matplotlib.pyplot import loglog,xlabel,ylabel,hist,figure,savefig,semilogx, semilogy,plot, errorbar, rc, tight_layout
from numpy import log10, array, zeros, isnan, logical_and, nan, linspace
from scipy.stats import linregress, spearmanr
from scipy.optimize import curve_fit
from pandas import read_excel, read_csv, DataFrame

from olm.USGS.loadWaterQualityData import loadSiteListData
from olm.general import CtoK
from olm.calcite import calc_K_H,concCaEqFromPCO2, pwpRateTheory, pwpRatePascal, pwp_to_mm_yr, pwpRateFranci, solutionFromCaPCO2, pwpFromSolution

rc('font', size=18)
class_colors = {0:'black', 1:'blue', 2:'red', 3:'green',4:'cyan'}

def check_plots_dir(sitesDir):
    plotsDir = os.path.join(sitesDir, 'plots/')
    if not os.path.exists(plotsDir):
        try:
            os.makedirs(plotsDir)
        except os.error:
            sys.exit("Problem creating plot directory: " + plotsDir)

def get_good_indicies(seriesList):
    test=0.
    for series in seriesList:
        test = test+series#Addition will create Nans if any of the series contains a NaN
    return ~test.isnull()

def calc_site_pwp(sitedf, sitephreeqc, returnPCO2=False):
    have_data = False
    if ('Stream flow, mean. daily' in sitedf.columns) and ( 'Temperature, water' in sitedf.columns) and ('CO2_Molality' in sitephreeqc.columns) and ('Ca+2_Activity'  in sitephreeqc.columns) and ('CO2_Activity'  in sitephreeqc.columns) and ('H+_Activity'  in sitephreeqc.columns) and ('HCO3-_Activity'  in sitephreeqc.columns):
       have_data = True
    if have_data:
        #create dataframe subset
        subdf = DataFrame({
                'Q':sitedf['Stream flow, mean. daily'],
                'T_C':sitedf['Temperature, water'],
                'CO2':sitephreeqc.CO2_Molality,
                'a_Ca':sitephreeqc['Ca+2_Activity'],
                'a_H2CO3s':sitephreeqc['CO2_Activity'],
                'a_H':sitephreeqc['H+_Activity'],
                'a_HCO3':sitephreeqc['HCO3-_Activity']
                })
        #Clear out NaN values
        subdf = subdf.dropna()
        if subdf.size>0:
            #Average any duplicate indicies
            g = subdf.groupby(level=0)#group by duplicate indicies
            subdf = g.mean()#average duplicate indicies
            #Calculate PCO2
            T_K = CtoK(subdf.T_C)
            K_H = calc_K_H(T_K)
            PCO2 = subdf.CO2/K_H    
            pwp_rates = pwpRateTheory(a_Ca=subdf.a_Ca, a_H2CO3s=subdf.a_H2CO3s, a_H=subdf.a_H, a_HCO3=subdf.a_HCO3, T_K=T_K, PCO2=PCO2)
            if returnPCO2:
                return [pwp_rates, PCO2]
            else:
                return pwp_rates
        else:
            if returnPCO2:
                return [None,None]
            else: 
                return None
    #Data were missing
    else:
        if returnPCO2:
            return [None,None]
        else: 
            return None
def pwp_Q_fit(Q, k, a, beta):
    F = k*(1. - a*Q**(-beta))
    return F

def make_pwp_vs_Q_plots(sitesDir, siteList=None, siteFile=None, makeHistograms=True, classFile='', bracket_charge_balance=False,plotloglog=False, nbins=10):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    if classFile != '':
        class_xls = read_excel(classFile, 'Sheet1',index_col=1, names=['name','site', 'recharge', 'age'])
    site_class = []
    rs =[]
    ps = []
    spp = []
    spr = []
    slopes= []
    plotsDir = os.path.join(sitesDir, 'plots/')
    #Check for plots directory and make it if one doesn't exist
    check_plots_dir(sitesDir)
    for site in list(sitePanelDict.keys()):
        if classFile !='':
            #read in classification for this site
            this_class = class_xls['recharge'][site]
            this_color=class_colors[this_class]
            site_class.append(this_class)
        else:
            this_color='black'
        print(("Making plot for: "+site))
        Q = sitePanelDict[site].data['Stream flow, mean. daily']
        pwp_rates = calc_site_pwp(sitePanelDict[site], sitePhreeqcDict[site])
        pwp_rates = pwp_to_mm_yr(pwp_rates)
        df = DataFrame({'Q':Q, 'pwp':pwp_rates})
        df = df.dropna()
        Qposids = df.Q>0
        df = df[Qposids]
        #Calculate fitting function parameters
        fit_converged = True
        if df.Q.values.size>5:
            try:
                opt_params, pcov = curve_fit(pwp_Q_fit, df['Q'], df['pwp'], p0=[1.,1.,1.])
                k = opt_params[0]
                a = opt_params[1]
                beta = opt_params[2]
            except RuntimeError as e:
                fit_converged = False
                print(e)
        else:
            fit_converged = False
        #Make Discharge vs. Saturation Ratio Plot
        figure()
        if plotloglog:
            if bracket_charge_balance:
                loglog(df.Q, df.pwp, '.', ms=0, color=this_color)
            else:
                loglog(df.Q, df.pwp, '.', color=this_color)
        else:
            if bracket_charge_balance:
                semilogx(df.Q, df.pwp, '.', color=this_color, ms=0)
            else:
                semilogx(df.Q, df.pwp, '.', color=this_color)
        #plot fitting function
        if fit_converged:
            Q_fit = linspace(df.Q.min(), df.Q.max(), 1000)
            pwp_fit = pwp_Q_fit(Q_fit, k, a, beta)
            plot(Q_fit, pwp_fit, '-k')
        xlabel('Discharge (cfs)')
        ylabel('PWP Dissolution Rate (mm/yr)')
        if bracket_charge_balance:
            #Load charge bracket data
            df_alk = read_csv(os.path.join(sitesDir, site, site+'-PHREEQC-Alk.csv'), parse_dates = True, index_col=0)
            df_ca = read_csv(os.path.join(sitesDir, site, site+'-PHREEQC-Ca.csv'), parse_dates = True, index_col=0)
            pwp_alk = pwp_to_mm_yr(calc_site_pwp(sitePanelDict[site], df_alk))
            pwp_ca = pwp_to_mm_yr(calc_site_pwp(sitePanelDict[site], df_ca))
            delta_alk = df.pwp - pwp_alk
            delta_ca = df.pwp - pwp_ca
            delta_upper = zeros(delta_ca.size)
            delta_lower = zeros(delta_ca.size)
            for i, delta in enumerate(delta_ca):
                if delta<0:
                    delta_upper[i] = abs(delta_ca[i])
                    delta_lower[i] = abs(delta_alk[i])
                elif delta>0:
                    delta_lower[i] = abs(delta_ca[i])
                    delta_upper[i] = abs(delta_alk[i])
            df['pwp_upper'] = delta_upper[Qposids.values]
            df['pwp_lower'] = delta_lower[Qposids.values]
            df['Alkalinity_balance_error'] = df_alk['Percent error']
            not_converged = df['Alkalinity_balance_error']>5.0
            df.pwp[not_converged]=nan
            df = df.dropna()
            df = df[df.Q>0] #USGS database has some negative values? -999999
            errorbar(df.Q, df.pwp, yerr=[df.pwp_upper, df.pwp_lower], fmt='o', color=this_color)

        tight_layout()
        savefig(os.path.join(plotsDir,site+'-Q_vs_PWP.pdf'))
        #Calculate correlation coefficient and regression
        slope,intercept,r,p,stderr = linregress(log10(df.Q), df.pwp)        
        rs.append(r)
        if isnan(r):
            return df
        ps.append(p)
        [this_spr, this_spp] = spearmanr(log10(df.Q), df.pwp)
        spr.append(this_spr)
        spp.append(this_spp)
        slopes.append(slope)
        if makeHistograms:
            figure()
            hist(df.pwp, normed=True, color=this_color)
            xlabel('PWP Rate mm/yr')
            ylabel('Frequency')
            savefig(os.path.join(plotsDir,site+'-PWPHist.pdf'))

    #Make histogram of pearson r values    
    figure()
    arr_rs = array(rs)
    arr_ps = array(ps)
    arr_spr = array(spr)
    arr_spp = array(spp)
    arr_site_class = array(site_class)
    goodids = logical_and(~isnan(arr_rs), ~isnan(arr_ps))
    arr_rs = arr_rs[goodids]
    arr_ps = arr_ps[goodids]
    arr_site_class = arr_site_class[goodids]
    if arr_rs[arr_site_class==1].shape[0]>1:
        arr1_rs_sig = arr_rs[logical_and(arr_site_class==1, arr_ps<0.05)]
        print("Class 1: number of sig ps = ",arr1_rs_sig.size,  "number insig = ", arr_rs[arr_site_class==1].size - arr1_rs_sig.size)       
        hist(arr1_rs_sig, nbins, color=class_colors[1], alpha=0.5)
    if arr_rs[arr_site_class==2].shape[0]>1:
        arr2_rs_sig = arr_rs[logical_and(arr_site_class==2, arr_ps<0.05)]
        print("Class 2: number of sig ps = ",arr2_rs_sig.size,  "number insig = ", arr_rs[arr_site_class==2].size - arr2_rs_sig.size)       
        hist(arr2_rs_sig, nbins, color=class_colors[2], alpha=0.5)
    if arr_rs[arr_site_class==3].shape[0]>1:
        arr3_rs_sig = arr_rs[logical_and(arr_site_class==3, arr_ps<0.05)]
        print("Class 3: number of sig ps = ",arr3_rs_sig.size,  "number insig = ", arr_rs[arr_site_class==3].size - arr3_rs_sig.size)       
        hist(arr3_rs_sig, nbins, color=class_colors[3], alpha=0.5)
    xlabel('Pearson R')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir,'PearsonR_Q_PWP.pdf'))
    #Make histogram of pearson p values
    figure()
    if arr_ps[arr_site_class==1].shape[0]>1:
        hist(arr_ps[arr_site_class==1], nbins, color=class_colors[1], alpha=0.5)
    if arr_ps[arr_site_class==2].shape[0]>1:
        hist(arr_ps[arr_site_class==2], nbins, color=class_colors[2], alpha=0.5)
    if arr_ps[arr_site_class==3].shape[0]>1:
        hist(arr_ps[arr_site_class==3], nbins, color=class_colors[3], alpha=0.5)
    xlabel('Pearson p')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir, 'PearsonP_Q_PWP.pdf'))

    #Make histogram of spearman r values
    arr_spr = array(spr)
    arr_spp = array(spp)
    arr_site_class = array(site_class)
    goodids = logical_and(~isnan(arr_spr), ~isnan(arr_spp))
    arr_spr = arr_spr[goodids]
    arr_spp = arr_spp[goodids]
    arr_site_class = arr_site_class[goodids]
    figure()
    if arr_spr[arr_site_class==1].shape[0]>1:
        arr1_spr_sig = arr_spr[logical_and(arr_site_class==1, arr_spp<0.05)]
        print("Class 1: number of sig ps = ",arr1_spr_sig.size,  "number insig = ", arr_spr[arr_site_class==1].size - arr1_spr_sig.size)       
        hist(arr1_spr_sig, nbins, color=class_colors[1], alpha=0.5)
    if arr_spr[arr_site_class==2].shape[0]>1:
        arr2_spr_sig = arr_spr[logical_and(arr_site_class==2, arr_spp<0.05)]
        print("Class 2: number of sig ps = ",arr2_spr_sig.size,  "number insig = ", arr_spr[arr_site_class==2].size - arr2_spr_sig.size)       
        hist(arr2_spr_sig, nbins, color=class_colors[2], alpha=0.5)
    if arr_spr[arr_site_class==3].shape[0]>1:
        arr3_spr_sig = arr_spr[logical_and(arr_site_class==3, arr_spp<0.05)]
        print("Class 3: number of sig ps = ",arr3_spr_sig.size,  "number insig = ", arr_spr[arr_site_class==3].size - arr3_spr_sig.size)       
        hist(arr3_spr_sig, nbins, color=class_colors[3], alpha=0.5)
    xlabel('Spearman R')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir, 'SpearmanR_Q_PWP.pdf'))
    #Make histogram of spearman p values
    figure()
    if arr_spp[arr_site_class==1].shape[0]>1:
        hist(arr_spp[arr_site_class==1], nbins, color=class_colors[1], alpha=0.5)
    if arr_spp[arr_site_class==2].shape[0]>1:
        hist(arr_spp[arr_site_class==2], nbins, color=class_colors[2], alpha=0.5)
    if arr_spp[arr_site_class==3].shape[0]>1:
        hist(arr_spp[arr_site_class==3], nbins, color=class_colors[3], alpha=0.5)
    xlabel('Spearman p')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir, 'SpearmanP_Q_PWP.pdf'))

def make_pwp_vs_saturation_ratio_plots(sitesDir, siteList=None, siteFile=None,  classFile='', plotloglog=True):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    if classFile != '':
        class_xls = read_excel(classFile, 'Sheet1',index_col=1, names=['name','site', 'recharge', 'age'])
    site_class = []
    rs =[]
    ps = []
    spp = []
    spr = []
    slopes= []
    plotsDir = os.path.join(sitesDir, 'plots/')
    #Check for plots directory and make it if one doesn't exist
    check_plots_dir(sitesDir)
    for site in list(sitePanelDict.keys()):
        if classFile !='':
            #read in classification for this site
            this_class = class_xls['recharge'][site]
            this_color=class_colors[this_class]
            site_class.append(this_class)
        else:
            this_color='black'
        print(("Making plot for: "+site))
        Q = sitePanelDict[site].data['Stream flow, mean. daily']
        pwp_rates = calc_site_pwp(sitePanelDict[site], sitePhreeqcDict[site])
        pwp_rates = pwp_to_mm_yr(pwp_rates)
        T_C = sitePanelDict[site].data['Temperature, water']
        T_K = CtoK(T_C)
        K_H = calc_K_H(T_K)
        CO2 = sitePhreeqcDict[site].CO2_Molality
        PCO2 = CO2/K_H
        CaEq = concCaEqFromPCO2(PCO2, T_C=T_C)
        sat_ratio = sitePhreeqcDict[site].Ca/CaEq
        df = DataFrame({'Q':Q, 'pwp':pwp_rates, 'sat_ratio':sat_ratio})
        df = df.dropna()
        #Make PWP vs. Saturation Ratio Plot
        figure()
        if plotloglog:
            loglog(df.pwp, df.sat_ratio, '.', color=this_color)
        else:
            semilogx(df.pwp, df.sat_ratio, '.', color=this_color)
        xlabel('PWP Dissolution Rate (mm/yr)')
        ylabel('Saturation ratio')
        savefig(os.path.join(plotsDir,site+'-Sat_ratio_vs_PWP.pdf'))
        #Calculate correlation coefficient and regression
        slope,intercept,r,p,stderr = linregress(df.pwp, df.sat_ratio)        
        rs.append(r)
        ps.append(p)
        [this_spr, this_spp] = spearmanr(df.pwp, df.sat_ratio)
        spr.append(this_spr)
        spp.append(this_spp)
        slopes.append(slope)

    #Make histogram of pearson r values    
    figure()
    arr_rs = array(rs)
    print("size arr_rs=", arr_rs.size)
    arr_ps = array(ps)
    arr_spr = array(spr)
    arr_spp = array(spp)
    arr_site_class = array(site_class)
    hist(arr_rs[arr_site_class==1],color=class_colors[1], alpha=0.5)
    hist(arr_rs[arr_site_class==2],color=class_colors[2], alpha=0.5)
    hist(arr_rs[arr_site_class==3],color=class_colors[3], alpha=0.5)
    xlabel('Pearson R')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir,'PearsonR_PWP_SatRatio.pdf'))
    #Make histogram of pearson p values
    figure()
    hist(arr_ps[arr_site_class==1],color=class_colors[1], alpha=0.5)
    hist(arr_ps[arr_site_class==2],color=class_colors[2], alpha=0.5)
    hist(arr_ps[arr_site_class==3],color=class_colors[3], alpha=0.5)
    xlabel('Pearson p')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir, 'PearsonP_PWP_SatRatio.pdf'))
    #Make histogram of spearman r values
    figure()
    hist(arr_spr[arr_site_class==1],color=class_colors[1], alpha=0.5)
    hist(arr_spr[arr_site_class==2],color=class_colors[2], alpha=0.5)
    hist(arr_spr[arr_site_class==3],color=class_colors[3], alpha=0.5)
    xlabel('Spearman R')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir, 'SpearmanR_PWP_SatRatio.pdf'))
    #Make histogram of spearman p values
    figure()
    hist(arr_spp[arr_site_class==1],color=class_colors[1], alpha=0.5)
    hist(arr_spp[arr_site_class==2],color=class_colors[2], alpha=0.5)
    hist(arr_spp[arr_site_class==3],color=class_colors[3], alpha=0.5)
    xlabel('Spearman p')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir, 'SpearmanP_PWP_SatRatio.pdf'))



def make_pwp_vs_T_plots(sitesDir, siteList=None, siteFile=None, makeHistograms=True, classFile=''):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    if classFile != '':
        class_xls = read_excel(classFile, 'Sheet1',index_col=1, names=['name','site', 'recharge', 'age'])
    rs = []
    ps = []
    slopes = []
    plotsDir = os.path.join(sitesDir, 'plots/')
    for site in list(sitePanelDict.keys()):
        if classFile !='':
            #read in classification for this site
            this_class = class_xls['recharge'][site]
            this_color=class_colors[this_class]
        else:
            this_color='black'
        print(("Making plot for: "+site))
        #Read out data into shorter variable names
        Q = sitePanelDict[site].data['Stream flow, mean. daily']
        T_C = sitePanelDict[site].data['Temperature, water']
        CO2 = sitePhreeqcDict[site].CO2_Molality
        Ca = sitePhreeqcDict[site]['Ca']
        a_Ca = sitePhreeqcDict[site]['Ca+2_Activity']
        a_H2CO3s = sitePhreeqcDict[site]['CO2_Activity']
        a_HCO3 = sitePhreeqcDict[site]['HCO3-_Activity']
        a_H = sitePhreeqcDict[site]['H+_Activity']
        #Filter out nan values
        #good_values = ~(Q+T_C+CO2+Ca).isnull()#Addition will create Nans if any of the series contains a NaN
        good_values = get_good_indicies([Q,T_C,CO2,Ca,a_Ca,a_H2CO3s,a_H, a_HCO3])
        Q = Q[good_values]
        T_C = T_C[good_values]
        CO2 = CO2[good_values]
        Ca = Ca[good_values]
        a_Ca = a_Ca[good_values]
        a_H = a_H[good_values]
        a_H2CO3s = a_H2CO3s[good_values]
        a_HCO3 = a_HCO3[good_values]
        #Calculate PCO2
        T_K = CtoK(T_C)
        K_H = calc_K_H(T_K)
        PCO2 = CO2/K_H
        #Calculate PWP rates
        pwp_rates = []
        for i, value in enumerate(Q):
#            sol = solutionFromCaPCO2(Ca[i], PCO2[i], T_C = T_C[i])
#            this_pwp = pwpFromSolution(sol)
            this_pwp =  pwpRateTheory(a_Ca=a_Ca[i], a_H2CO3s=a_H2CO3s[i], a_H=a_H[i], a_HCO3=a_HCO3[i], T_K=T_K[i], PCO2=PCO2[i])
            this_rate_mm_yr = pwp_to_mm_yr(this_pwp)
            pwp_rates.append(this_rate_mm_yr)
        #Check for plots directory and make it if one doesn't exist
        check_plots_dir(sitesDir)
        #Make Discharge vs. Saturation Ratio Plot
        figure()
        plot(T_C, pwp_rates, 'o', color=this_color)
        xlabel('Temperature ($^\circ C$)')
        ylabel('PWP Dissolution Rate (mm/yr)')
        savefig(os.path.join(plotsDir,site+'-T_vs_PWP.pdf'))
        #Calculate correlation coefficient and regression
        slope,intercept,r,p,stderr = linregress(T_C, pwp_rates)
        rs.append(r)
        ps.append(p)
        slopes.append(slope)
        if makeHistograms:
            figure()
            hist(pwp_rates, normed=True, color=this_color)
            xlabel('PWP Rate mm/yr')
            ylabel('Frequency')
            savefig(os.path.join(plotsDir,site+'-PWPHist.pdf'))

    #Make histogram of pearson r values    
    figure()
    hist(rs,normed=True)
    xlabel('Pearson R')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir,'PearsonR_T_PWP.pdf'))


def make_saturation_index_vs_Q_plots(sitesDir, siteList=None, siteFile=None, makeHistograms=True):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    rs = []
    ps = []
    slopes = []
    plotsDir = os.path.join(sitesDir, 'plots/')
    for site in list(sitePanelDict.keys()):
        print(("Making plot for: "+site))
        #Read out data into shorter variable names
        Q = sitePanelDict[site].data['Stream flow, mean. daily']
        SI = sitePhreeqcDict[site].SI_Calcite
        #Filter out nan values
        good_values = get_good_indicies([Q,SI])
        Q = Q[good_values]
        SI = SI[good_values]
        check_plots_dir(sitesDir)
        #Make Discharge vs. Saturation Index Plot
        figure()
        semilogx(Q, SI, 'o')
        xlabel('Discharge (cfs)')
        ylabel('$SI_{calcite}$')
        savefig(os.path.join(plotsDir,site+'-Q_vs_SI.pdf'))
        #Calculate correlation coefficient and regression
        slope,intercept,r,p,stderr = linregress(log10(Q), SI)
        rs.append(r)
        ps.append(p)
        slopes.append(slope)
        if makeHistograms:
            figure()
            hist(SI, normed=True)
            xlabel('$SI_{calcite}$')
            ylabel('Frequency')
            savefig(os.path.join(plotsDir,site+'-SIHist.pdf'))        

    #Make histogram of pearson r values    
    figure()
    hist(rs,normed=True)
    xlabel('Pearson R')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir,'PearsonR_Q_SIRatio.pdf'))

def make_saturation_ratio_vs_Q_plots(sitesDir, siteList=None, siteFile=None, makeHistograms=True):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    rs = []
    ps = []
    slopes = []
    plotsDir = os.path.join(sitesDir, 'plots/')
    for site in list(sitePanelDict.keys()):
        print(("Making plot for: "+site))
        #Read out data into shorter variable names
        Q = sitePanelDict[site].data['Stream flow, mean. daily']
        T_C = sitePanelDict[site].data['Temperature, water']
        CO2 = sitePhreeqcDict[site].CO2_Molality
        Ca = sitePhreeqcDict[site]['Ca']
        #Filter out nan values
        #good_values = ~(Q+T_C+CO2+Ca).isnull()#Addition will create Nans if any of the series contains a NaN
        good_values = get_good_indicies([Q,T_C,CO2,Ca])
        Q = Q[good_values]
        T_C = T_C[good_values]
        CO2 = CO2[good_values]
        Ca = Ca[good_values]
        #Calculate PCO2
        T_K = CtoK(T_C)
        K_H = calc_K_H(T_K)
        PCO2 = CO2/K_H
        #Calculate Equilibrium Ca concentration
        CaEq = concCaEqFromPCO2(PCO2, T_C=T_C)
        #Check for plots directory and make it if one doesn't exist
        check_plots_dir(sitesDir)
        #Make Discharge vs. Saturation Ratio Plot
        figure()
        loglog(Q, Ca/CaEq, 'o')
        xlabel('Discharge (cfs)')
        ylabel('Saturation Ratio $[Ca]/[Ca]_{eq}$')
        savefig(os.path.join(plotsDir,site+'-Q_vs_SatRatio.pdf'))
        #Calculate correlation coefficient and regression
        slope,intercept,r,p,stderr = linregress(log10(Q), Ca/CaEq)
        rs.append(r)
        ps.append(p)
        slopes.append(slope)
        if makeHistograms:
            figure()
            hist(Ca/CaEq, normed=True)
            xlabel('Saturation Ratio $[Ca]/[Ca]_{eq}$')
            ylabel('Frequency')
            savefig(os.path.join(plotsDir,site+'-SatRatioHist.pdf'))

    #Make histogram of pearson r values    
    figure()
    hist(rs,normed=True)
    xlabel('Pearson R')
    ylabel('Frequency')
    savefig(os.path.join(plotsDir,'PearsonR_Q_SatRatio.pdf'))



def make_pwp_vs_PCO2_plots(sitesDir, siteList=None, siteFile=None, makeHistograms=True, classFile='', bracket_charge_balance=False,plotloglog=False, nbins=10):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    if classFile != '':
        class_xls = read_excel(classFile, 'Sheet1',index_col=1, names=['name','site', 'recharge', 'age'])
    site_class = []
    rs =[]
    ps = []
    spp = []
    spr = []
    slopes= []
    plotsDir = os.path.join(sitesDir, 'plots/')
    #Check for plots directory and make it if one doesn't exist
    check_plots_dir(sitesDir)
    for site in list(sitePanelDict.keys()):
        if classFile !='':
            #read in classification for this site
            this_class = class_xls['recharge'][site]
            this_color=class_colors[this_class]
            site_class.append(this_class)
        else:
            this_color='black'
        print(("Making plot for: "+site))
        Q = sitePanelDict[site].data['Stream flow, mean. daily']
        pwp_rates,PCO2 = calc_site_pwp(sitePanelDict[site], sitePhreeqcDict[site], returnPCO2=True)
        pwp_rates = pwp_to_mm_yr(pwp_rates)
        df = DataFrame({'PCO2':PCO2, 'pwp':pwp_rates})
        df = df.dropna()
        #Make PWP vs. PCO2 Plot
        figure()
        if plotloglog:
            if bracket_charge_balance:
                loglog(df.PCO2, df.pwp, '.', ms=0, color=this_color)
            else:
                loglog(df.PCO2, df.pwp, '.', color=this_color)
        else:
            if bracket_charge_balance:
                semilogx(df.PCO2, df.pwp, '.', color=this_color, ms=0)
            else:
                semilogx(df.PCO2, df.pwp, '.', color=this_color)
        xlabel('PCO2 (atm)')
        ylabel('PWP Dissolution Rate (mm/yr)')
        if bracket_charge_balance:
            #Load charge bracket data
            df_alk = read_csv(os.path.join(sitesDir, site, site+'-PHREEQC-Alk.csv'), parse_dates = True, index_col=0)
            df_ca = read_csv(os.path.join(sitesDir, site, site+'-PHREEQC-Ca.csv'), parse_dates = True, index_col=0)
            pwp_alk = pwp_to_mm_yr(calc_site_pwp(sitePanelDict[site], df_alk))
            pwp_ca = pwp_to_mm_yr(calc_site_pwp(sitePanelDict[site], df_ca))
            delta_alk = df.pwp - pwp_alk
            delta_ca = df.pwp - pwp_ca
            delta_upper = zeros(delta_ca.size)
            delta_lower = zeros(delta_ca.size)
            for i, delta in enumerate(delta_ca):
                if delta<0:
                    delta_upper[i] = abs(delta_ca[i])
                    delta_lower[i] = abs(delta_alk[i])
                elif delta>0:
                    delta_lower[i] = abs(delta_ca[i])
                    delta_upper[i] = abs(delta_alk[i])
            df['pwp_upper'] = delta_upper
            df['pwp_lower'] = delta_lower
            df['Alkalinity_balance_error'] = df_alk['Percent error']
            not_converged = df['Alkalinity_balance_error']>5.0
            df.pwp[not_converged]=nan
            df = df.dropna()
            errorbar(df.PCO2, df.pwp, yerr=[df.pwp_upper, df.pwp_lower], fmt='o', color=this_color)

        savefig(os.path.join(plotsDir,site+'-PCO2_vs_PWP.pdf'))
        #Calculate correlation coefficient and regression
#        slope,intercept,r,p,stderr = linregress(log10(df.Q), df.pwp)        
#        rs.append(r)
#        if isnan(r):
#            return df
#        ps.append(p)
#        [this_spr, this_spp] = spearmanr(log10(df.Q), df.pwp)
#        spr.append(this_spr)
#        spp.append(this_spp)
#        slopes.append(slope)
#        if makeHistograms:
#            figure()
#            hist(df.pwp, normed=True, color=this_color)
#            xlabel('PWP Rate mm/yr')
#            ylabel('Frequency')
#            savefig(os.path.join(plotsDir,site+'-PWPHist.pdf'))

    #Make histogram of pearson r values    
#    figure()
#    arr_rs = array(rs)
#    arr_ps = array(ps)
#    arr_spr = array(spr)
#    arr_spp = array(spp)
#    arr_site_class = array(site_class)
#    goodids = logical_and(~isnan(arr_rs), ~isnan(arr_ps))
#    arr_rs = arr_rs[goodids]
#    arr_ps = arr_ps[goodids]
#    arr_site_class = arr_site_class[goodids]
#    if arr_rs[arr_site_class==1].shape[0]>1:
#        arr1_rs_sig = arr_rs[logical_and(arr_site_class==1, arr_ps<0.05)]
#        print "Class 1: number of sig ps = ",arr1_rs_sig.size,  "number insig = ", arr_rs[arr_site_class==1].size - arr1_rs_sig.size       
#        hist(arr1_rs_sig, nbins, color=class_colors[1], alpha=0.5)
#    if arr_rs[arr_site_class==2].shape[0]>1:
#        arr2_rs_sig = arr_rs[logical_and(arr_site_class==2, arr_ps<0.05)]
#        print "Class 2: number of sig ps = ",arr2_rs_sig.size,  "number insig = ", arr_rs[arr_site_class==2].size - arr2_rs_sig.size       
#        hist(arr2_rs_sig, nbins, color=class_colors[2], alpha=0.5)
#    if arr_rs[arr_site_class==3].shape[0]>1:
#        arr3_rs_sig = arr_rs[logical_and(arr_site_class==3, arr_ps<0.05)]
#        print "Class 3: number of sig ps = ",arr3_rs_sig.size,  "number insig = ", arr_rs[arr_site_class==3].size - arr3_rs_sig.size       
#        hist(arr3_rs_sig, nbins, color=class_colors[3], alpha=0.5)
#    xlabel('Pearson R')
#    ylabel('Frequency')
#    savefig(os.path.join(plotsDir,'PearsonR_Q_PWP.pdf'))
#    #Make histogram of pearson p values
#    figure()
#    if arr_ps[arr_site_class==1].shape[0]>1:
#        hist(arr_ps[arr_site_class==1], nbins, color=class_colors[1], alpha=0.5)
#    if arr_ps[arr_site_class==2].shape[0]>1:
#        hist(arr_ps[arr_site_class==2], nbins, color=class_colors[2], alpha=0.5)
#    if arr_ps[arr_site_class==3].shape[0]>1:
#        hist(arr_ps[arr_site_class==3], nbins, color=class_colors[3], alpha=0.5)
#    xlabel('Pearson p')
#    ylabel('Frequency')
#    savefig(os.path.join(plotsDir, 'PearsonP_Q_PWP.pdf'))

#    #Make histogram of spearman r values
#    arr_spr = array(spr)
#    arr_spp = array(spp)
#    arr_site_class = array(site_class)
#    goodids = logical_and(~isnan(arr_spr), ~isnan(arr_spp))
#    arr_spr = arr_spr[goodids]
#    arr_spp = arr_spp[goodids]
#    arr_site_class = arr_site_class[goodids]
#    figure()
#    if arr_spr[arr_site_class==1].shape[0]>1:
#        arr1_spr_sig = arr_spr[logical_and(arr_site_class==1, arr_spp<0.05)]
#        print "Class 1: number of sig ps = ",arr1_spr_sig.size,  "number insig = ", arr_spr[arr_site_class==1].size - arr1_spr_sig.size       
#        hist(arr1_spr_sig, nbins, color=class_colors[1], alpha=0.5)
#    if arr_spr[arr_site_class==2].shape[0]>1:
#        arr2_spr_sig = arr_spr[logical_and(arr_site_class==2, arr_spp<0.05)]
#        print "Class 2: number of sig ps = ",arr2_spr_sig.size,  "number insig = ", arr_spr[arr_site_class==2].size - arr2_spr_sig.size       
#        hist(arr2_spr_sig, nbins, color=class_colors[2], alpha=0.5)
#    if arr_spr[arr_site_class==3].shape[0]>1:
#        arr3_spr_sig = arr_spr[logical_and(arr_site_class==3, arr_spp<0.05)]
#        print "Class 3: number of sig ps = ",arr3_spr_sig.size,  "number insig = ", arr_spr[arr_site_class==3].size - arr3_spr_sig.size       
#        hist(arr3_spr_sig, nbins, color=class_colors[3], alpha=0.5)
#    xlabel('Spearman R')
#    ylabel('Frequency')
#    savefig(os.path.join(plotsDir, 'SpearmanR_Q_PWP.pdf'))
#    #Make histogram of spearman p values
#    figure()
#    if arr_spp[arr_site_class==1].shape[0]>1:
#        hist(arr_spp[arr_site_class==1], nbins, color=class_colors[1], alpha=0.5)
#    if arr_spp[arr_site_class==2].shape[0]>1:
#        hist(arr_spp[arr_site_class==2], nbins, color=class_colors[2], alpha=0.5)
#    if arr_spp[arr_site_class==3].shape[0]>1:
#        hist(arr_spp[arr_site_class==3], nbins, color=class_colors[3], alpha=0.5)
#    xlabel('Spearman p')
#    ylabel('Frequency')
#    savefig(os.path.join(plotsDir, 'SpearmanP_Q_PWP.pdf'))


def make_PCO2_vs_Q_plots(sitesDir, siteList=None, siteFile=None, classFile='', makeHistograms=True):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    if classFile != '':
        class_xls = read_excel(classFile, 'Sheet1',index_col=1, names=['name','site', 'recharge', 'age'])
    rs = []
    ps = []
    site_class = []
    slopes = []
    plotsDir = os.path.join(sitesDir, 'plots/')
    for site in list(sitePanelDict.keys()):
        if classFile !='':
            #read in classification for this site
            this_class = class_xls['recharge'][site]
            this_color=class_colors[this_class]
            site_class.append(this_class)
        else:
            this_color='black'
        print(("Making plot for: "+site))
        #Read out data into shorter variable names
        Q = sitePanelDict[site].data['Stream flow, mean. daily']
        T_C = sitePanelDict[site].data['Temperature, water']
        CO2 = sitePhreeqcDict[site].CO2_Molality
        Ca = sitePhreeqcDict[site]['Ca']
        #Filter out nan values
        #good_values = ~(Q+T_C+CO2+Ca).isnull()#Addition will create Nans if any of the series contains a NaN
        good_values = get_good_indicies([Q,T_C,CO2,Ca])
        Q = Q[good_values]
        T_C = T_C[good_values]
        CO2 = CO2[good_values]
        Ca = Ca[good_values]
        #Calculate PCO2
        T_K = CtoK(T_C)
        K_H = calc_K_H(T_K)
        PCO2 = CO2/K_H
        figure()
        loglog(Q, PCO2, 'o', color=this_color)
        xlabel('Discharge (cfs)')
        ylabel('PCO2 (atm)')
        savefig(os.path.join(plotsDir,site+'-PCO2_vs_Q.pdf'))
        #Calculate correlation coefficient and regression
#        slope,intercept,r,p,stderr = linregress(log10(Q), Ca/CaEq)
#        rs.append(r)
#        ps.append(p)
#        slopes.append(slope)
#        if makeHistograms:
#            figure()
#            hist(Ca/CaEq, normed=True)
#            xlabel('Saturation Ratio $[Ca]/[Ca]_{eq}$')
#            ylabel('Frequency')
#            savefig(os.path.join(plotsDir,site+'-SatRatioHist.pdf'))

#    #Make histogram of pearson r values    
#    figure()
#    hist(rs,normed=True)
#    xlabel('Pearson R')
#    ylabel('Frequency')
#    savefig(os.path.join(plotsDir,'PearsonR_Q_SatRatio.pdf'))


def make_PCO2_vs_T_plots(sitesDir, siteList=None, siteFile=None, classFile='', makeHistograms=True):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    if classFile != '':
        class_xls = read_excel(classFile, 'Sheet1',index_col=1, names=['name','site', 'recharge', 'age'])
    site_class = []
    plotsDir = os.path.join(sitesDir, 'plots/')
    for site in list(sitePanelDict.keys()):
        if classFile !='':
            #read in classification for this site
            this_class = class_xls['recharge'][site]
            this_color=class_colors[this_class]
            site_class.append(this_class)
        else:
            this_color='black'
        print(("Making plot for: "+site))
        #Read out data into shorter variable names
        T_C = sitePanelDict[site].data['Temperature, water']
        CO2 = sitePhreeqcDict[site].CO2_Molality
        good_values = get_good_indicies([T_C,CO2])
        T_C = T_C[good_values]
        CO2 = CO2[good_values]
        #Calculate PCO2
        T_K = CtoK(T_C)
        K_H = calc_K_H(T_K)
        PCO2 = CO2/K_H
        figure()
        loglog(T_C, PCO2, 'o', color=this_color)
        xlabel('Temperature (C)')
        ylabel('PCO2 (atm)')
        savefig(os.path.join(plotsDir,site+'-PCO2_vs_T.pdf'))

def make_PCO2_vs_Month_plots(sitesDir, siteList=None, siteFile=None, classFile=''):
    sitePanelDict, sitePhreeqcDict = loadSiteListData(processedSitesDir=sitesDir, loadPhreeqc=True)
    if classFile != '':
        class_xls = read_excel(classFile, 'Sheet1',index_col=1, names=['name','site', 'recharge', 'age'])
    site_class = []
    plotsDir = os.path.join(sitesDir, 'plots/')
    for site in list(sitePanelDict.keys()):
        if classFile !='':
            #read in classification for this site
            this_class = class_xls['recharge'][site]
            this_color=class_colors[this_class]
            site_class.append(this_class)
        else:
            this_color='black'
        print(("Making plot for: "+site))
        #Read out data into shorter variable names
        T_C = sitePanelDict[site].data['Temperature, water']
        CO2 = sitePhreeqcDict[site].CO2_Molality
        good_values = get_good_indicies([T_C,CO2])
        T_C = T_C[good_values]
        CO2 = CO2[good_values]
        #Calculate PCO2
        T_K = CtoK(T_C)
        K_H = calc_K_H(T_K)
        PCO2 = CO2/K_H
        sample_month = PCO2.index.month
        figure()
        #semilogy(sample_month, PCO2, 'o', color=this_color)
        plot(sample_month, PCO2, 'o', color=this_color)
        xlabel('Month')
        ylabel('PCO2 (atm)')
        savefig(os.path.join(plotsDir,site+'-PCO2_vs_Month.pdf'))
