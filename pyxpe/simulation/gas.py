#!/usr/bin/env python
# Copyright (C) 2007--2016 the X-ray Polarimetry Explorer (XPE) team.
#
# For the license terms see the file LICENSE, distributed along with this
# software.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Main class for running a simulation job (an `experiment')
# follows roughly the concept of xpesim/TExperiment

import os
import numpy as np
from scipy.interpolate import interp1d
from pyxpe.logging_ import logger 

# get dir name of this module
import pyxpe.simulation
THIS_MODULE_DIR = os.path.dirname(pyxpe.simulation.__file__)

MIX_FILE_NAME   = os.path.join(THIS_MODULE_DIR,'gasproperties/MIXTURES.DAT')
COMP_FILE_NAME  = os.path.join(THIS_MODULE_DIR,'gasproperties/COMPOUNDS.DAT')
ELE_FILE_NAME   = os.path.join(THIS_MODULE_DIR,'gasproperties/ELEMENTS.DAT')
XSEC_FILE_NAME  = os.path.join(THIS_MODULE_DIR,'gasproperties/xsections/')

#0.082057338[L atm K-1 mol-1]*(273.15 + 25)[K]/1[atm] https://en.wikipedia.org/wiki/Molar_volume
MOLAR_VOLUME    = 24465.395 # at 25 C and 1 atm in cm #used to be 22414. valid at 0C
AVOGADRO_NUMBER = 6.02e23
MIN_PHOTOELECTRIC_DATA_ENERGY      = 1.0
MAX_PHOTOELECTRIC_DATA_ENERGY      = 1000.0
MIN_ANALYTIC_CROSS_SECTIONS_ENERGY = 0.05

class gasmix:
    """Gas mixture object
    """
    def __init__(self, mixid, pressure, mix_file_name = MIX_FILE_NAME):
        self.mixid    = mixid
        self.pressure = pressure
        logger.info("Mix Id = %d Pressure %s (atm)" % (mixid, pressure))
        self.parse_mix_file(mix_file_name)
        self.eval_properties()
        
    def parse_mix_file(self, mix_file_name):
        """ Parsing Mixture file and get relevant properties
        """
        mix_file = file(mix_file_name)
        mix_line = mix_file.readlines()[self.mixid+1].split()
        self.comp_list = [compound(mix_line[0], float(mix_line[1]))]
        self.diffusion_sigma = float(mix_line[6])
        self.diffusion_sigma /= np.sqrt(self.pressure) # IS THIS TRUE?
        if mix_line[2] != '-':
            self.comp_list.append(compound(mix_line[2], float(mix_line[3])))
        if mix_line[4] != '-':
            self.comp_list.append(compound(mix_line[4], float(mix_line[5])))
        logger.info("Mix of %d gas, Diffusion Sigma %.3f um (for 1 cm drift)"%\
                     (len(self.comp_list),self.diffusion_sigma))
        for comp in self.comp_list:
            logger.info(str(comp))

    def eval_properties(self):
        """
        """
        TotalStopPower = 0.
        TotalSecondary = 0.
        for comp in self.comp_list:
            TotalStopPower += comp.getfraction()*comp.get_spoppingpower()
            TotalSecondary += comp.getfraction()*comp.get_ionsnumber()
        # ionization
        self.WIonization = TotalStopPower/TotalSecondary;
        logger.info("W ionization %f" %self.WIonization)

        # create list of elements
        self.Density   = 0.0;
        self.elements_list  = [] # create list of elements
        self.elements_names = [] 
        for comp in self.comp_list:
            for j in xrange(comp.get_nelements()):
                myElement = element(comp.get_element(j))
                #Evaluate Mixture density
                mydensity = myElement.AtomicWeight*\
                            comp.get_atoms(j)*\
                            comp.getfraction()*\
                            self.pressure/MOLAR_VOLUME
                self.Density += mydensity
                ElemName = myElement.ChemicalSymbol
                if ElemName not in self.elements_names:
                    myElement.Density = mydensity
                    self.elements_list.append(myElement)
                    self.elements_names.append(ElemName)
                else:
                    elementID = self.elements_names.index(ElemName)
                    self.elements_list[elementID].Density += mydensity
                    
        logger.info("Density of the mixture= %.8f g/cm^3" % self.Density)
        for elem in self.elements_list:
            logger.info(str(elem))

    def GetElectronRange(self, energy):
        """\brief Returns the electron range as a function of the energy.
        Energy to be provided in keV.
        """
        if isinstance( energy , (int, float) ):
            energy = np.array([energy])
        return (10.0/self.Density)*np.power(10, (-5.1 + 1.358*np.log10(energy) + 0.215*(np.power(np.log10(x),2)) -0.043*(np.power(np.log10(x),3))))

    def GetPhotoelectricCrossSection(self, energy):
        """\brief Returns the Photoelectric cross section at a given energy. 
        Energy to be provided in keV, cross section returned in cm^2/g.
        """
        # min/max energy check TBD
        CrossSection = 0;
        for elem in self.elements_list:
            CrossSection +=elem.GetPhotoelectricCrossSection(energy)*elem.Density/self.Density
        return CrossSection
            
    def GetAbsorptionLength(self, energy):
        """\brief Returns the photoelectric absorption lenght at a given energy
        Energy to be provided in keV, absorption lenght returned in cm.
        """
        return 1.0/(self.Density*self.GetPhotoelectricCrossSection(energy))
    
    def GetEfficiency(self, energy, gap_tickness):
        """ \brief Returns the detection efficiency at given energy 
        and gap thickness.
        Energy to be provided in keV; gap tickness in mm
        """
        # notice that GetAbsorptionLength out is in cm
        return (1 -np.exp(-0.1*gap_tickness/self.GetAbsorptionLength(energy)))
    
    def GetConvertingElement(self, energy, r0):
        """ \brief EXTRACT the element converting an incoming photon.
        The element in the mixture is extracted according to 
        the relative values of the photoelectric cross sections 
        at a given energy.
        Energy to be provided in keV.
        """
        total_x_sect = self.GetPhotoelectricCrossSection(energy)*self.Density
        partial_probability = 0
        for elem in self.elements_list:
            absorption_prob = elem.GetPhotoelectricCrossSection(energy)*elem.Density/total_x_sect
            if (r0 < partial_probability+absorption_prob):
                return elem
            partial_probability += absorption_prob
        # this should never happen!
        logger.error("Converting element not extracted!!!")
        return None

        
    def GetElasticMeanFreePath(self, energy, mode):
        """\brief Returns the the total elestic mean free path for the mixture.
        The mean free path is evaluated by summing the inverse of the mean free
        paths for all the elements.
        Energy to be provided in keV, mode cam be either RUTHERFORD or MOTT.
        """
        if mode.lower() != "rutherford" and mode.lower() != "mott":
            logger.error("Mode can be only rutherford or mott, not %s" %\
                          mode)
            return -1
        MeanFreePath = 0.
        for elem in self.elements_list:
            MeanFreePath += 1./(elem.GetElasticMeanFreePath(energy,
                                                            elem.Density,
                                                            mode))
        return 1./MeanFreePath
        
    def GetScatteringElement(self, energy, mode, r0):
        """ \brief EXTRACTthe element scattering a travelling photoelectron.
        The element in the mixture is extracted according to the relative 
        values ofthe elastic cross sections at a given energy.
        Energy to be provided in keV, mode can be either RUTHERFORD OR MOTT.
        """
        total_x_sect = self.GetElasticMeanFreePath(energy, mode) # 1/x-section
        partial_probability = 0
        for elem in self.elements_list:
            scat_prob = total_x_sect/elem.GetElasticMeanFreePath(energy,
                                                                 elem.Density,
                                                                 mode)
            if (r0 < partial_probability+scat_prob):
                return elem
            partial_probability += scat_prob
        # this should never happen!
        logger.error("Scattering element not extracted!!!")
        return None
        
    def GetStoppingPower(self, energy):
        """\brief Returns the stopping power at a given energy.
        Energy to be provided in keV, outcome in keV/cm.
        """
        stopping_power = 0.
        for elem in self.elements_list:
            stopping_power += elem.GetStoppingPower(energy)
        return stopping_power


"""
CLASS ELEMENT
"""
class element:
    def __init__(self, name, ele_file_name = ELE_FILE_NAME):
        ele_file  = file(ele_file_name)
        ele_lines = ele_file.readlines()
        self.ChemicalSymbol = None
        for ll in ele_lines:
            l = ll.split()
            if len(l)>0 and l[0] == name:
                self.ChemicalSymbol    = name
                self.AtomicNumber      =  int(l[1])
                self.AtomicWeight      =  float(l[2])
                self.kEdge             =  float(l[3])
                self.FluorescenceYield =  float(l[4])
                break
        if self.ChemicalSymbol == None:
            logger.error("No element %s found, please update file %s" %\
                          (name, ele_file_name))

        # The Density is re-calculated in gasmix:eval_properties
        # as partial density
        self.Density = self.AtomicWeight/MOLAR_VOLUME;  
  
        # Returns the mean ionization potential.
        # This parametrization is from Berger and Seltzer (1964);
        # see Joy's book for the reference.
        # The ionization potential is returned in keV.
        self.MeanIonizationPotential = 0.001*(9.76*self.AtomicNumber +\
                                              58.5/np.power(self.AtomicNumber, 0.19))
        # eval X section
        self.GetPhotoelectricCrossSection = self.EvalPhotoelectricCrossSection()
        
        
    def __str__(self):
        return "Element %s Z=%d A= %.2f partial density=%.8f g/cm^3" %\
            (self.ChemicalSymbol, self.AtomicNumber, self.AtomicWeight,
             self.Density)

    
    def EvalPhotoelectricCrossSection(self):
        """\brief Evaluates the photoelectric cross section. 
        Reads the files containing the photoelectric cross section for a given
        energy grid (one per element, taken from the NIST database).
        Energy to be provided in keV, cross section returned in cm^2/g. 
        """
        xsec_file  = file(os.path.join(XSEC_FILE_NAME,
                                       self.ChemicalSymbol+'.txt'))
        xsec_lines = xsec_file.readlines()
        energies  = []
        PhotoelectricCrossSection = []
        # assume this marker always exist
        beginId = xsec_lines.index('BEGIN:\r\n')+1
        xsec_lines = xsec_lines[beginId:]
        for l in xsec_lines:
            ll = l.split()
            if len(ll)>0:
                energies.append(float(ll[0])*1000.0) # keV
                PhotoelectricCrossSection.append(float(ll[3]))

        self.energies = np.array(energies)
        self.PhotoelectricCrossSection = np.array(PhotoelectricCrossSection)
        # return and interpolate 1D object, diff with spline interp is 1%
        return interp1d(np.array(energies), np.array(PhotoelectricCrossSection))
                    
    def GetPhotoelectricCrossSection(self, energy):
        return self.PhotoelectricCrossSection(energy)
    
    def GetRutherfordScreeningFactor(self, energy):
        """ \brief Returns the Rutherford screening factor at a given anergy.
        This parametrization is from Bishop (1976); see Joy's book for the
        reference.                   
        Energy to be provided in keV (accurate down to 50 eV)
        """
        if isinstance( energy , (int, float) ):
            energy = np.array([energy])
        if energy.any()< MIN_ANALYTIC_CROSS_SECTIONS_ENERGY:
            logger.error("Cannot eval screening factor below %f keV"%\
                          MIN_ANALYTIC_CROSS_SECTIONS_ENERGY)
            return -1
        return (3.4E-3)*np.power(self.AtomicNumber, 0.67)/energy
    
    def GetRutherfordTotalCrossSection(self, energy):
        """\brief Returns the Rutherford integral cross section.  
        This parametrization is given to Newbury and Myklebust (1981); see
        Joy's book for the reference.  
        Energy to be provided in keV, cross section returned in cm^2/atom 
        (accurate down to 50 eV). 
        """
        if isinstance( energy , (int, float) ):
            energy = np.array([energy])
        if energy.any()< MIN_ANALYTIC_CROSS_SECTIONS_ENERGY:
            logger.error("Cannot eval Rutherford X section below %f keV"%\
                          MIN_ANALYTIC_CROSS_SECTIONS_ENERGY)
            return -1
        
        alpha = self.GetRutherfordScreeningFactor(energy);
        xsection = 5.21e-21
        xsection *= (np.power(self.AtomicNumber, 2.0)/np.power(energy, 2.0))
        xsection *= (4*np.pi/(alpha*(1 + alpha)))
        xsection *= np.power(((energy + 511)/(energy + 1024)), 2.0)
        return xsection
    
    def GetMottTotalCrossSection(self, energy):
        """\brief Returns the Mott integral cross section at a given energy.
        This parametrization is given to Browning (1992); see Joy's book 
        for the reference
        Energy to be provided in keV, cross section returned in cm^2/atom.
        (accurate down to 50 eV).
        """
        if isinstance( energy , (int, float) ):
            energy = np.array([energy])
        if energy.any()< MIN_ANALYTIC_CROSS_SECTIONS_ENERGY:
            logger.error("Cannot eval Mott X section below %f keV"%\
                          MIN_ANALYTIC_CROSS_SECTIONS_ENERGY)
            return -1
        
        u = np.log10(8*energy*np.power(self.AtomicNumber, -1.33))
        out = (4.7E-18)*(np.power(self.AtomicNumber, 1.33) + 0.032*np.power(self.AtomicNumber, 2.0))/((energy + 0.0155*np.power(self.AtomicNumber, 1.33)*np.power(energy, 0.5))*(1 - 0.02*np.power(self.AtomicNumber, 0.5)*np.exp(-u*u)))
        return out
    
    def GetElasticMeanFreePath(self, energy, density, mode):
        """\brief Returns the mean free path for elastic collisions 
        at a given energy.
        Energy to be provided in keV, density in g/cm^3; mode can be either
        RUTHERFORD or MOTT. The mean free path is returned in cm.
        """
        if isinstance( energy , (int, float) ):
            energy = np.array([energy])
        if energy.any()< MIN_ANALYTIC_CROSS_SECTIONS_ENERGY:
            logger.error("Cannot eval ElasticMeanFreePath below %f keV"%\
                          MIN_ANALYTIC_CROSS_SECTIONS_ENERGY)
            return -1

        MeanFreePath = self.AtomicWeight/(AVOGADRO_NUMBER*density)
        if mode.lower()=='rutherford':
            MeanFreePath /= self.GetRutherfordTotalCrossSection(energy);
            return MeanFreePath;
        elif mode.lower()=='mott':
            MeanFreePath /= self.GetMottTotalCrossSection(energy);
            return MeanFreePath;
        else:
            logger.error("Mode can be only rutherford or mott, not %s" %\
                          mode)
            return -1
        
    def GetScatteringAngle(self, energy, mode, r0, r1):
        """Returns a RANDOM scattering angle at a given energy.
        Energy to be provided in keV; mode can be either RUTHERFORD or MOTT.
        The cosine of the angle is returned.
        r0 and r1 are 2 random number (from external generator)
        """
        if isinstance( energy , (int, float) ):
            energy = np.array([energy])
        if energy.any()< MIN_ANALYTIC_CROSS_SECTIONS_ENERGY:
            logger.error("Cannot eval Scattering Angle below %f keV"%\
                          MIN_ANALYTIC_CROSS_SECTIONS_ENERGY)
            return -1
        screening_factor = self.GetRutherfordScreeningFactor(energy)
        if mode.lower()=='rutherford':
            return np.arccos(1.-(2.*screening_factor*r0)/(1.+screening_factor-r1))
        elif mode.lower()=='mott':
            CorrectionFactor = 2.2 - ((92.0 - self.AtomicNumber)/92.0);
            return np.arccos(1.-(2.*screening_factor*np.power((r0), CorrectionFactor))/(1.+screening_factor-r1))
        else:
            logger.error("Mode can be only rutherford or mott, not %s" %\
                          mode)
            return -1
        
    def GetStoppingPower(self, energy):
        """ \brief Returns the stopping power at a given energy.
        This parametrization is from Joy and Luo (1989); 
        see Joy's book for the reference.
        Energy to be provided in keV, outcome in keV/cm.
        """
        if isinstance( energy , (int, float) ):
            energy = np.array([energy])
        if energy.any()< MIN_ANALYTIC_CROSS_SECTIONS_ENERGY:
            logger.error("Cannot eval StoppingPower below %f keV"%\
                          MIN_ANALYTIC_CROSS_SECTIONS_ENERGY)
            return -1
        stop_pwr = 78500*self.AtomicNumber*np.log(1.166*(energy + 0.85*self.MeanIonizationPotential)/self.MeanIonizationPotential)/(energy*self.AtomicWeight)*self.Density
        return stop_pwr
        
"""
CLASS COMPOUND aka molecule with fractional part
"""     
class compound:
    def __init__(self, name, fraction, comp_file_name = COMP_FILE_NAME):
        """ Init compound and eval properties from table
        """
        self.name = name
        self.fraction = fraction
        comp_file = file(comp_file_name)
        comp_lines = comp_file.readlines()
        self.n_elements = 0
        self.elements = []
        self.atoms    = []        
        for ll in comp_lines:
            l = ll.split()
            if l[0] == self.name:
                self.n_elements = int(l[1])
                self.NIonsT     = float(l[2]) # total ionization/cm
                self.StoppingP  = float(l[3]) # de/dx in ev/cm
                for i in range(4, len(l), 2):
                    self.elements.append(l[i])
                    self.atoms.append(int(l[i+1]))
        if self.n_elements == 0:
            logger.error("No compound %s found, please update file %s" %\
                          (self.name, comp_file_name))

    def getname(self):
        return self.name
    
    def getfraction(self):
        return self.fraction
    
    def get_spoppingpower(self):
        """GetCompoundStoppingPower
        dE/dx ev/cm
        """
        return self.StoppingP

    def get_ionsnumber(self):
        """GetCompoundIonsNumber
        total ions/cm
        """
        return self.NIonsT

    def get_nelements(self):
        """GetnElementsInCompound
        """
        return self.n_elements

    def get_element(self, i):
        """GetElementsInCompound
        Element name of the element i
        """
        return self.elements[i]

    def get_atoms(self, i):
        """GetnAtomsInCompound 
        Num of atoms of the element i
        """
        return self.atoms[i]
        
    def __str__(self):
        return "Compound %s fraction %f" % (self.name, self.fraction)
        



    

if __name__ == '__main__':
    
    print "------------------------"
    g0 = gasmix(26, 0.8)
    print "------------------------"
    g = gasmix(12, 1.0)
    #print "------------------------"
    #g1 = gasmix(16, 1.0)
    
    
    import matplotlib.pyplot as plt
    xnew = np.logspace(0,3, num=101)

    #plt.plot(xnew, g.GetPhotoelectricCrossSection(xnew), '.', label='PhXsec')
    #plt.plot(xnew, g0.GetPhotoelectricCrossSection(xnew), '.', label='PhXsec0')
    #plt.plot(xnew, g.GetAbsorptionLength(xnew), '.', label='absLen')
    #plt.plot(xnew, g0.GetAbsorptionLength(xnew), '.', label='absLen0')
    #plt.plot(xnew, g.GetEfficiency(xnew, 10), '.', label='eff10')
    #plt.plot(xnew, g.GetEfficiency(xnew, 100), '.', label='eff100')
    
    #plt.xlabel('E (keV)')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.legend()
    #
    plt.show()
        
    el = []
    for ename in ['He', 'C','O','H']:
        el.append(element(ename))
        for energy in [3.7, 5.9]:
            print ename, energy, "-", \
                el[-1].GetPhotoelectricCrossSection(energy),\
                el[-1].GetRutherfordTotalCrossSection(energy),\
                el[-1].GetMottTotalCrossSection(energy),\
                el[-1].GetElasticMeanFreePath(energy, 1, "mott")
        
        plt.plot(xnew, el[-1].GetPhotoelectricCrossSection(xnew), '.',\
                 label = "%s"%ename)
        """
        plt.plot(xnew, el[-1].GetRutherfordTotalCrossSection(xnew), '.',\
                 label = "Ruth")
        plt.plot(xnew, el[-1].GetMottTotalCrossSection(xnew), '.',\
                 label = "Mott")
        """
    #plt.ylabel(ename)
    plt.ylabel('Photoelectrinc x-section')
    plt.xlabel('E (keV)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
        
        

    """
    import matplotlib.pyplot as plt
    el = []
    for ename in ['He', 'C','O','H']:
        el.append(element(ename))
        plt.plot(el[-1].energies, el[-1].PhotoelectricCrossSection, 'o', label = ename)

    plt.xlabel('E (keV)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
    
    """

    """
    e = element('C')
    print e
    # test interpolation method
    x = e.energies
    y = e.PhotoelectricCrossSection
    xnew = np.logspace(0,3, num=201)
    
    import matplotlib.pyplot as plt
    plt.plot(x, y, 'o', label='data')
    
    from scipy.interpolate import interp1d
    f = interp1d(x, y)
    plt.plot(xnew, f(xnew), '-', label='interp1d')
    #f1 = interp1d(x, y, kind='cubic') # una merda!
    #plt.plot(xnew, f1(xnew), '-', label='interp1d cubic')

    from scipy import interpolate
    tck = interpolate.splrep(x, y, s=0)
    ynew = interpolate.splev(xnew, tck, der=0)
    plt.plot(xnew, ynew, '-', label='spline cubic')

    ff1 = plt.figure(1)
    plt.xlabel('E (keV)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='best')

    ff2 = plt.figure(2)
    plt.xlabel('E (keV)')
    plt.xscale('log')
    plt.plot(xnew, ynew/f(xnew), '.')
    plt.show()
    """
