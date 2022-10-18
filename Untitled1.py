#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyopenms


# In[3]:


# Avogadro number
print("Avogadro's number is: ",pyopenms.Constants.AVOGADRO)


# In[4]:


#ElementDB
from pyopenms import *
e=ElementDB()
e.hasElement('H')


# In[57]:


e.hasElement('v')


# In[28]:


#Nitrogen element
nitrogen=e.getElement('N')
print("Name is: ",nitrogen.getName())
print("MonoWeight of Nirtrogen is: ",nitrogen.getMonoWeight())
print("Nitrogen's symbol is: ",nitrogen.getSymbol())
print("Nitrogen's AverageWeight is: ",nitrogen.getAverageWeight())

print("one mole of nitrogen weights: ",2*nitrogen.getAverageWeight(),"grams")
isotopes=nitrogen.getIsotopeDistribution()
print(isotopes) # print as an object


# In[32]:


#isotope
#Nirogen
print("Nitrogen: ")
nitrogen_isoDist={"mass": [], "abundance": []}
for i in isotopes.getContainer():
    print("Nitrogen isotope: ",i.getMZ(),",Has abundance: ",i.getIntensity()*100,"%")
    nitrogen_isoDist["mass"].append(i.getMZ())
    nitrogen_isoDist["abundance"].append((i.getIntensity()*100))
    
    


print("Hydrogen: ")    
#Hydrogen  
hydrogen_isoDist = {"mass": [],"abundance": []}
hydrogen=e.getElement('H')
isotopes=hydrogen.getIsotopeDistribution()

for i in isotopes.getContainer():
    print("hydrogen isotope: ",i.getMZ(),",Has abundance: ",i.getIntensity()*100,"%")
    hydrogen_isoDist["mass"].append(i.getMZ())
    hydrogen_isoDist["abundance"].append((i.getIntensity()*100))
    


# In[33]:


#display isotopes distribution of Nitrogen
import math
from matplotlib import pyplot as plt

def adjustText(x1, y1, x2, y2):
    if y1 > y2:
        plt.annotate('%0.3f' % (y2), xy=(x2, y2), xytext=(x2+0.5,y2+9),
                     textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
                     horizontalalignment='right', verticalalignment='top')
    else:
        plt.annotate('%0.3f' % (y1), xy=(x1, y1), xytext=(x1+0.5,y1+9),
                     textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
                     horizontalalignment='right', verticalalignment='top')
        
def plotDistribution(distribution):
    n = len(distribution["mass"])
    for i in range(0, n):
        plt.vlines(x=distribution["mass"][i], ymin=0, ymax=distribution["abundance"][i])
        if int(distribution["mass"][i - 1]) == int(distribution["mass"][i])                 and i != 0:
            adjustText(distribution["mass"][i - 1], distribution["abundance"][i - 1],
                       distribution["mass"][i], distribution["abundance"][i])
        else:
            plt.text(x=distribution["mass"][i],
                     y=(distribution["abundance"][i] + 2),
                     s='%0.3f' % (distribution["abundance"][i]), va='center',
                     ha='center')
    plt.ylim([0, 110])
    plt.xticks(range(math.ceil(distribution["mass"][0]) - 2,
                     math.ceil(distribution["mass"][-1]) + 2))

#Nitrogen figure   
plt.figure(figsize=(12,8)) 
plt.subplot(1,2,1)
plt.title("Isotopic distribution of nitrogen")
plotDistribution(nitrogen_isoDist)
plt.xlabel("Atomic mass (am)")
plt.ylabel("Relative abundance (%)")

#Hydrogen figure
plt.figure(figsize=(12,8))
plt.subplot(1,2,2)
plt.title("Isotopic distribution of hydrogen")
plotDistribution(hydrogen_isoDist)
plt.xlabel("Atomic mass (am)")
plt.ylabel("Relative abundance (%)")
plt.show()


# In[37]:


#Mass Defect
oxygen_isotopes=e.getElement('O').getIsotopeDistribution().getContainer()

oxygen_isotope_difference=oxygen_isotopes[1].getMZ() - oxygen_isotopes[0].getMZ()


hydrogen_isotope=e.getElement('H').getIsotopeDistribution().getContainer()

hydrogen_isotope_difference=hydrogen_isotope[1].getMZ()-hydrogen_isotope[0].getMZ()


print ("Mass difference between 17O and 16O:", oxygen_isotope_difference)

print ("Mass difference between 15H and 14H:", hydrogen_isotope_difference)

print("Relative deviation: ", 100*(hydrogen_isotope_difference - oxygen_isotope_difference)/hydrogen_isotope_difference ,"%")


# In[39]:


#Molecular Formula
methanol=EmpiricalFormula("CH3OH")

water=EmpiricalFormula("H2O")

ethanol=EmpiricalFormula("CH2")+methanol
#print ethanol formula
print("Ethanol chemical formula: ",ethanol.toString())
#print ethanol composition
print("Ethanol Composition from:",ethanol.getElementalComposition())

#print number of specific element in ethanol 
print("Ethanol has: ",ethanol.getElementalComposition()[b'H'], "hydrogen atoms")


# In[40]:


#isotopes

ethanol = EmpiricalFormula("C2H6O")
print("Ethanol chemical formula:", ethanol.toString())
print("Ethanol composition:", ethanol.getElementalComposition())
print("Ethanol weight:", ethanol.getMonoWeight())

ethanol = EmpiricalFormula("(13)C1CH6O")
print("Ethanol chemical formula:", ethanol.toString())
print("Ethanol composition:", ethanol.getElementalComposition())
print("Ethanol weight:", ethanol.getMonoWeight())

ethanol = EmpiricalFormula("(13)C2H6O")
print("Ethanol chemical formula:", ethanol.toString())
print("Ethanol composition:", ethanol.getElementalComposition())
print("Ethanol weight:", ethanol.getMonoWeight())

    


# In[43]:


#Amino Acids

pro=ResidueDB().getResidue("Proline")

#display name
print(pro.getName())

#display three letter code
print(pro.getThreeLetterCode())

#display one Letter code
print(pro.getOneLetterCode())

#display Average Weight
print(pro.getAverageWeight())

#display Mono weight
print(pro.getMonoWeight())

#display the value of Pka 
print(pro.getPka())

#display Formula
print(pro.getFormula().toString())


# In[50]:


#Amino Acid Modifications

oxid=ModificationsDB().getModification("Oxidation")
print(oxid.getUniModAccession())
print(oxid.getUniModRecordId())
print(oxid.getDiffMonoMass())
print(oxid.getId())
print(oxid.getFullId())
print(oxid.getFullName())
print(oxid.getDiffFormula())


# In[52]:


#investigate the isotopic distribution of the modification 
#the isotopic pattern of the modification (Oxygen)
isotopes=oxid.getDiffFormula().getIsotopeDistribution(CoarseIsotopePatternGenerator(5))
for i in isotopes.getContainer():
    print(i.getMZ(), " : ",i.getIntensity())


# In[56]:


#Ribonucleotides
cysteine=RibonucleotideDB().getRibonucleotide(b"C")
print(cysteine.getName())
print(cysteine.getCode())
print(cysteine.getAvgMass())
print(cysteine.getMonoMass())
print(cysteine.getFormula().toString())
print(cysteine.isModified())
methyladenosine = RibonucleotideDB().getRibonucleotide(b"m6A")
print(methyladenosine.getName())
print(methyladenosine.isModified())


# In[ ]:




