# -*- coding: utf-8 -*-
"""
Silvija Milosavljevic
"""
### individual-based model of polyploid establishment in a population 
### code used for Milosavljevic et al. 2024 American Journal of Botany
### A metabolic perspective on polyploid invasion and the emergence of  
### life histories: insights from a mechanistic model

### spatial model - diploids and tetraploids in a population, 
### depend on size to consume nutrients,
### reproductive energy-dependent number of offspring, 
### offspring is randomly dispersed and that is only movement,
### reproduction happens after a season, only once in life,
### prob of offspring becoming tetraploid is sampled from a beta distr, 
### mu is the average polyploid formation 
### code includes visuals for initial inspection of the simulation results

import random as rnd
import numpy as np
import matplotlib as mtl
import matplotlib.pyplot as plt
from scipy.stats import beta

class Individual:
    '''Class that regulates individuals and their properties'''
    def __init__(self,
                 x,
                 y,
                 ploidy, 
                 base_weight, # weight at birth
                 birthday, 
                 parent, 
                 Energy_Reserve):
        '''Initialization'''
        self.x = x
        self.y = y
        self.ploidy = ploidy
        self.current_weight = base_weight # weight at birth
        self.Energy_Reserve = Energy_Reserve # what is left over after somatic maintenance and growth
        self.birthday = birthday
        self.parent = parent # parent ploidy so that it is trackable, it can be 2 or 4

class Metapopulation:
    '''Contains the whole population, regulates daily affairs'''
    def __init__(self, 
                 max_x, 
                 max_y):
        '''Initialization'''           
        self.max_x = max_x
        self.max_y = max_y
        initial_resources = 10 # initial and maximum for each position in the space 
        self.environment = np.zeros((self.max_x,self.max_y)) + initial_resources
        self.population = []
        self.initialize_pop()

        self.x=[]  #timer
        self.y1=[] #diploids
        self.y2=[] #tetraploids
        self.y3=[] #resources
        self.y4=[] #proportion of tetraploids

        self.cwd = [] # current weight of diploids 
        self.cwt = [] # current weight of tetraploids
        self.erd = [] # energy for reproduction of diploids
        self.ert = [] # energy for reproduction of tetraploids

        self.tetrafromdi = [] # number of newly born 4x from 2x parents
        self.tetrafromtetra = [] # number of newly born 4x from 4x parents

        self.maturityDi = []
        self.maturityTetra = [] # for tracking the average maturity time for each timestep

        self.weightatreprDi = []
        self.weightatreprTetra = [] # for tracking the weight at reproduction

    def initialize_pop(self):
        '''Initialize individuals'''
        startpop = 1600 
        # pop size at the start is 1600, 1 individual per cell on average
        # for each starting individual, random spot in space is chosen
        # weight is sampled from random distribution for seeds
        for n in range(startpop):
            x = rnd.randint(0,self.max_x-1)
            y = rnd.randint(0,self.max_y-1)
            baseweight = np.random.normal(0.5,0.05) 
            # all are seeds in the start
            self.population.append(Individual(x, y, 2, baseweight, 0, 2, 0)) 
            # all inds are 2x at start

    def a_day_in_the_life(self, timer):
        
        # empty distributions for plotting diploids and tetraploids in space
        Dist_pop_di=np.zeros((self.max_x, self.max_y))
        Dist_pop_tetra=np.zeros((self.max_x, self.max_y))            
        rnd.shuffle(self.population)
        # shuffle population so that individuals in the beginning of the list
        # don't get an advantage 
        # at the end of all the processes, counting of diploids and tetraploids 
        # for the tracking and figures
        countDiploids, countTetraploids, fromdiploid, fromtetraploid = 0, 0, 0, 0
        for indiv in self.population:
            if indiv.ploidy == 2:
                countDiploids += 1
            elif indiv.ploidy == 4:
                countTetraploids += 1
                if indiv.parent == 2:
                    fromdiploid += 1
                elif indiv.parent == 4:
                    fromtetraploid += 1
        current_w_di = []
        current_w_tetra = []
        current_er_di = []
        current_er_tetra = []
        timetomaturityDi = []
        timetomaturityTetra = []
        weightatreproductionDi = []
        weightatreproductionTetra = []
        # printing during the simulation itself: 
        print(countDiploids, countTetraploids, file = output_file)
        print(fromdiploid, fromtetraploid, file= output_file2)

        oldpop = self.population[:]
        print('Gen:  '+str(timer)+'  Population: '+str(len(oldpop)))
        del self.population[:]

        for indiv in oldpop:

            #if timer % 100 == 0: 
                # every 100 generations write out all the indivs and their properties, 
                # or every generation if needed: 
            output_file1.write('%20s'%str(timer) + '\t'+ \
                    '%20s'%str(indiv.x) + '\t'+ \
                        '%20s'%str(indiv.y) + '\t'+ \
                            '%20s'%str(indiv.current_weight) + '\t'+ \
                                '%20s'%str(indiv.Energy_Reserve) + '\t'+ \
                                    '%20s'%str(indiv.ploidy) + '\t'+ \
                                        '%20s'%str(indiv.birthday) + '\t'+\
                                            '%20s'%str(indiv.parent) + '\t'+\
                                                '%20s'%str(self.environment[int(indiv.x), int(indiv.y)]) + '\n') 
            if indiv.ploidy == 2: 
                current_w_di.append(indiv.current_weight)
                current_er_di.append(indiv.Energy_Reserve)
            elif indiv.ploidy == 4: 
                current_w_tetra.append(indiv.current_weight)
                current_er_tetra.append(indiv.Energy_Reserve)

            background_mortality = rnd.uniform(0,1)

            if  timer == 0 or timer % 100 != 0:

                if background_mortality > 0.0005: # 5 in 10k individuals dies due to chance

                    # Photosynthesis, resource partitioning 
                    if indiv.current_weight > 0:
                        PhotosynRate = (32.8 * (indiv.current_weight+indiv.Energy_Reserve) / (1.06 + 0.94 * (indiv.current_weight+indiv.Energy_Reserve)) ) * 60 * 60 * 24 * 10**(-6)
                        
                        if indiv.ploidy == 2: 
                            BMR = 0.36
                            VegGrowthR = 0.54
                            ReprGrowthR = 0.1
                        elif indiv.ploidy == 4:
                            BMR = 0.361
                            VegGrowthR = 0.54
                            ReprGrowthR = 0.099
                            # FOR DIPLOIDS: 
                                # 10% of photosynthesis is allocated to reproduction
                                # 54% of photosynthesis is allocated to somatic weight growth
                                # 36% is lost to basal metabolism
                            # FOR TETRAPLOIDS: 
                                # CHANGE BMR, VegGrowthR AND ReprGrowthR ACCORDINGLY
                        if self.environment[int(indiv.x), int(indiv.y)] > 0:
                            if self.environment[int(indiv.x), int(indiv.y)] >= PhotosynRate:
                                self.environment[int(indiv.x), int(indiv.y)] -= PhotosynRate
                                indiv.Energy_Reserve += ReprGrowthR * PhotosynRate
                                indiv.current_weight += VegGrowthR * PhotosynRate

                            elif self.environment[int(indiv.x), int(indiv.y)] >= (1-ReprGrowthR)*PhotosynRate:
                                indiv.current_weight += VegGrowthR * PhotosynRate
                                if self.environment[int(indiv.x), int(indiv.y)] - (1-ReprGrowthR)*PhotosynRate > 0: 
                                    indiv.Energy_Reserve += (self.environment[int(indiv.x), int(indiv.y)] - (1-ReprGrowthR)*PhotosynRate)
                                self.environment[int(indiv.x), int(indiv.y)] = 0
                            elif BMR*PhotosynRate < self.environment[int(indiv.x), int(indiv.y)] < (1-ReprGrowthR)*PhotosynRate:
                                indiv.current_weight += self.environment[int(indiv.x), int(indiv.y)] - BMR*PhotosynRate
                                self.environment[int(indiv.x), int(indiv.y)] = 0
                            elif self.environment[int(indiv.x), int(indiv.y)] < BMR*PhotosynRate:
                                p = self.environment[int(indiv.x), int(indiv.y)]
                                while p < BMR*PhotosynRate and indiv.Energy_Reserve > 0:
                                    p += 0.01
                                    indiv.Energy_Reserve -= 0.01
                                if p < BMR* PhotosynRate: 
                                    indiv.current_weight -= (BMR* PhotosynRate - p)
                                self.environment[int(indiv.x), int(indiv.y)] = 0
                            # if the energy reserve is negative, force it back to 0 and lose vegetative w:
                            if indiv.Energy_Reserve < 0:
                                indiv.current_weight -= abs(indiv.Energy_Reserve)
                                indiv.Energy_Reserve = 0
                        else: # if env has no resources at that moment
                            if indiv.Energy_Reserve >= BMR * PhotosynRate: 
                                indiv.Energy_Reserve -= BMR * PhotosynRate
                            elif 0 < indiv.Energy_Reserve < BMR * PhotosynRate:
                                indiv.current_weight = indiv.current_weight - (BMR * PhotosynRate - indiv.Energy_Reserve)
                                indiv.Energy_Reserve = 0
                            elif indiv.Energy_Reserve == 0:
                                indiv.current_weight = indiv.current_weight - BMR * PhotosynRate
                            elif indiv.Energy_Reserve < 0:
                                indiv.current_weight -= abs(indiv.Energy_Reserve)
                                indiv.Energy_Reserve = 0
                            # did not get any resources but lost the metabolic costs from reserve and weight
                        # DO NOT FORGET TO APPEND THE INDIVIDUALS BACK INTO THE POPULATION 
                        self.population.append(indiv)
                        if indiv.ploidy == 2:
                            Dist_pop_di[int(indiv.x), int(indiv.y)]+=1
                        else:
                            Dist_pop_tetra[int(indiv.x), int(indiv.y)]+=1
                    else: 
                        continue # die from losing all weight
                else:
                    continue # dies due to background mortality (chance)
            
            # Reproduction on 100th day, if there is any reproductive reserve
            elif timer > 0 and timer % 100 == 0: 
                if indiv.Energy_Reserve > 0: 
                    if indiv.ploidy == 2:
                        output_file5.write('%20s'%str(timer) + '\t'+ \
                                                '%20s'%str((timer - indiv.birthday)/100) + '\t'+\
                                                    '%20s'%str(indiv.current_weight) + '\n') 
                        while indiv.Energy_Reserve > 0:
                            # while the energy for repr is available, try to make an offs
                            probUnreducedGamete = beta.rvs(a=2, b=82, size=1)
                            
                            probabilityPolyploid = rnd.uniform(0,1)
                            output_file9.write('%20s'%str(timer) + '\t'+ \
                                        '%20s'%str(indiv.x) + '\t'+ \
                                            '%20s'%str(indiv.y) + '\t'+ \
                                                '%20s'%str(probUnreducedGamete)  + '\t' +\
                                                    '%20s'%str(probabilityPolyploid) +'\n') 

                            ploidyLevel = 2
                            baseweight = np.random.normal(0.5,0.05)
                            if probabilityPolyploid <= probUnreducedGamete:
                                ploidyLevel = 4
                                baseweight = np.random.normal(0.5,0.05)
                            if (indiv.Energy_Reserve - baseweight) >= 0:
                                x_coords = []
                                y_coords = []
                                for i in range(0,self.max_x):
                                    for j in range(0, self.max_y):
                                        if i >= indiv.x - 1 or i <= indiv.x + 1:
                                            x_coords.append(i)
                                        if j >= indiv.y - 1 or j <= indiv.y + 1:
                                            y_coords.append(j)
                                x_offspring = x_coords[rnd.randint(0, len(x_coords)-1)]
                                y_offspring = y_coords[rnd.randint(0, len(y_coords)-1)]

                                self.population.append(Individual(x_offspring,
                                                                y_offspring,
                                                                ploidyLevel, 
                                                                baseweight,
                                                                timer, 2, 0))
                                if ploidyLevel == 2: 
                                    Dist_pop_di[int(x_offspring), int(y_offspring)]+=1
                                else: 
                                    Dist_pop_tetra[int(x_offspring), int(y_offspring)]+=1
                                output_file7.write('%20s'%str(timer) + '\t'+ \
                                        '%20s'%str(x_offspring) + '\t'+ \
                                            '%20s'%str(y_offspring) + '\t'+ \
                                                '%20s'%str(baseweight) + '\t'+ \
                                                    '%20s'%str(ploidyLevel) + '\n') 

                                indiv.Energy_Reserve = indiv.Energy_Reserve - baseweight
                                timetomaturityDi.append((timer - indiv.birthday)/100)
                                weightatreproductionDi.append(indiv.current_weight)
                                if timer % 100 == 0: 
                                    # every 100 generations write out all successfully reproducing indivs
                                    output_file3.write('%20s'%str(timer) + '\t'+ \
                                        '%20s'%str(indiv.x) + '\t'+ \
                                            '%20s'%str(indiv.y) + '\t'+ \
                                                '%20s'%str(indiv.current_weight) + '\t'+ \
                                                    '%20s'%str(indiv.Energy_Reserve) + '\t'+ \
                                                        '%20s'%str(indiv.ploidy) + '\t'+ \
                                                            '%20s'%str(indiv.birthday) + '\t'+\
                                                                '%20s'%str(indiv.parent) + '\t'+\
                                                                    '%20s'%str(self.environment[int(indiv.x), int(indiv.y)]) + '\n') 

                            else:
                                break
                    elif indiv.ploidy == 4:
                        output_file6.write('%20s'%str(timer) + '\t'+ \
                                                '%20s'%str((timer - indiv.birthday)/100) + '\t'+\
                                                    '%20s'%str(indiv.current_weight) + '\n')                         
                        while indiv.Energy_Reserve > 0:
                            ploidyLevel = 4
                            baseweight = np.random.normal(0.5,0.05) 
                            if (indiv.Energy_Reserve - baseweight) >= 0: 
                                # if the parent has enough energy for this offspring, 
                                # otherwise it is done 
                                x_coords = []
                                y_coords = []
                                for i in range(0,self.max_x):
                                    for j in range(0, self.max_y):
                                        if i >= indiv.x - 1 or i <= indiv.x + 1:
                                            x_coords.append(i)
                                        if j >= indiv.y - 1 or j <= indiv.y + 1:
                                            y_coords.append(j)
                                x_offspring = x_coords[rnd.randint(0, len(x_coords)-1)]
                                y_offspring = y_coords[rnd.randint(0, len(y_coords)-1)]
                                
                                self.population.append(Individual(x_offspring,
                                                                y_offspring,
                                                                ploidyLevel, 
                                                                baseweight,
                                                                timer, 4, 0))
                                Dist_pop_tetra[int(x_offspring), int(y_offspring)]+=1
                                output_file8.write('%20s'%str(timer) + '\t'+ \
                                        '%20s'%str(x_offspring) + '\t'+ \
                                            '%20s'%str(y_offspring) + '\t'+ \
                                                '%20s'%str(baseweight) + '\t'+ \
                                                    '%20s'%str(ploidyLevel) + '\n') 
                                indiv.Energy_Reserve = indiv.Energy_Reserve - baseweight
                                timetomaturityTetra.append((timer - indiv.birthday)/100)
                                weightatreproductionTetra.append(indiv.current_weight)
                                if timer % 100 == 0:
                                    output_file3.write('%20s'%str(timer) + '\t'+ \
                                        '%20s'%str(indiv.x) + '\t'+ \
                                            '%20s'%str(indiv.y) + '\t'+ \
                                                '%20s'%str(indiv.current_weight) + '\t'+ \
                                                    '%20s'%str(indiv.Energy_Reserve) + '\t'+ \
                                                        '%20s'%str(indiv.ploidy) + '\t'+ \
                                                            '%20s'%str(indiv.birthday) + '\t'+\
                                                                '%20s'%str(indiv.parent) + '\t'+\
                                                                    '%20s'%str(self.environment[int(indiv.x), int(indiv.y)]) + '\n') 
                            else:
                                break
                    '''elif indiv.Energy_Reserve <= 0 and indiv.ploidy == 4:
                    self.population.append(indiv)
                    Dist_pop_tetra[int(indiv.x), int(indiv.y)]+=1'''
                    # only 4x are perennials, for that scenario, unquote the previous 3 lines
                else:
                    self.population.append(indiv)
                    if indiv.ploidy == 2:
                        Dist_pop_di[int(indiv.x), int(indiv.y)]+=1
                    else:
                        Dist_pop_tetra[int(indiv.x), int(indiv.y)]+=1
                        # all indivs are perennials
                        # for all annuals scenario, comment previous 6 lines                   

        for x in range(0,self.max_x):
            for y in range(0,self.max_y):
                self.environment[x,y] += 1 
                np.clip(self.environment, 0, 10, out = self.environment)
        # amount of resources has to stay between 0 and 10

        self.x.append(timer)
        self.y1.append(countDiploids)
        self.y2.append(countTetraploids)
        self.y3.append(np.ndarray.mean(self.environment))
        if (countTetraploids+countDiploids) != 0:
            self.y4.append(countTetraploids/(countTetraploids+countDiploids))
        else:
            self.y4.append(0)
        self.cwd.append(np.mean(current_w_di))
        self.cwt.append(np.mean(current_w_tetra))
        self.erd.append(np.mean(current_er_di))
        self.ert.append(np.mean(current_er_tetra))

        self.maturityDi.append(np.mean(timetomaturityDi))
        self.maturityTetra.append(np.mean(timetomaturityTetra))

        self.weightatreprDi.append(np.mean(weightatreproductionDi))
        self.weightatreprTetra.append(np.mean(weightatreproductionTetra))

        if (fromdiploid+fromtetraploid) != 0:
            self.tetrafromdi.append(fromdiploid/(fromdiploid+fromtetraploid))
            self.tetrafromtetra.append(fromtetraploid/(fromdiploid+fromtetraploid))
        elif (fromdiploid+fromtetraploid) == 0: 
            self.tetrafromdi.append(0)
            self.tetrafromtetra.append(0)
        if timer > 0 and timer % 100 == 0: 
            faxmaturityDih.hist(timetomaturityDi, histtype='step', stacked=True, fill=False)
            faxmaturityTetrah.hist(timetomaturityTetra, histtype='step', stacked=True, fill=False)

        ################################################
        # Create figure final timestep of the simulation
        if timer == 9999:
            print(np.mean(self.tetrafromdi), np.mean(self.tetrafromtetra))
            print(np.mean(self.tetrafromdi), np.mean(self.tetrafromtetra), np.mean(self.y4), file = output_file4)
            print(np.mean(self.y4))

            fax1.imshow(Dist_pop_di, animated=False, cmap='Reds', interpolation='none', origin="upper")
            fax11.imshow(Dist_pop_tetra, animated=False, cmap='Blues', interpolation='none', origin="upper")
            fax2.imshow(self.environment, animated=False,vmax=100, cmap='YlOrBr', interpolation='none', origin="upper")
            fax3.plot(self.x, self.y1,'r', animated=False)
            fax4.plot(self.x, self.y2,'b', animated=False)
            fax5.plot(self.x, self.y3,'chocolate', animated=False)

            faxcwDh.hist(current_w_di, color="r", fill=False)
            faxcwTh.hist(current_w_tetra, color="b", fill=False)
            faxcwD.plot(self.x, self.cwd,'r', animated=False)
            faxcwT.plot(self.x, self.cwt,'b', animated=False)

            faxerDh.hist(current_er_di, color="r", fill=False)
            faxerTh.hist(current_er_tetra, color="b", fill=False)
            faxerD.plot(self.x, self.erd,'r', animated=False)
            faxerT.plot(self.x, self.ert,'b', animated=False)

            faxProp.plot(self.x, self.y4, "b", animated=False)
        
            ax.bar(self.x, self.tetrafromdi, label='4x from 2x parent', color="firebrick")
            ax.bar(self.x, self.tetrafromtetra, bottom=self.tetrafromdi, label='4x from 4x parent', color=(0.03137254901960784, 0.18823529411764706, 0.4196078431372549, 1.0))
            ax.set_ylabel('Proportion')
            ax.set_title('Proportion of 4x from 2x and 4x as parents')
            ax.legend()

            faxmaturityDi.plot(self.x, self.maturityDi,'r.', animated=False)
            faxmaturityTetra.plot(self.x, self.maturityTetra,'b.', animated=False)

            faxwatreprDi.plot(self.x, self.weightatreprDi, "r.", animated=False)
            faxwatreprTetra.plot(self.x, self.weightatreprTetra, "b.", animated=False)
        

iter = 0
maxIter = 9 # How many simulations to be repeated

while (iter <= maxIter):
    meta = Metapopulation(40,40)

    fig1, ((fax1, fax11), (fax3, fax4), (fax2, fax5)) = plt.subplots(3, 2, figsize=(7,7), constrained_layout=True)
    fig2, ((faxcwD, faxcwDh), (faxcwT, faxcwTh)) = plt.subplots(2, 2, figsize=(7,7), constrained_layout=True)
    fig4, ((faxerD, faxerDh), (faxerT, faxerTh)) = plt.subplots(2, 2, figsize=(7,7), constrained_layout=True)
    fig3, (faxProp) = plt.subplots(1,1, figsize=(7,7), constrained_layout=True)
    fig5, ax = plt.subplots()

    fig6, (faxmaturityDi, faxmaturityTetra) = plt.subplots(1, 2, figsize=(8,8), constrained_layout=True)
    fig7, (faxmaturityDih, faxmaturityTetrah) = plt.subplots(1, 2, figsize=(7,7), constrained_layout=True)
    fig8, (faxwatreprDi, faxwatreprTetra) = plt.subplots(1, 2, figsize=(7,7), constrained_layout=True)

    faxwatreprDi.set_title('Diploid weight at reproduction')
    faxwatreprTetra.set_title('Tetraploid weight at reproduction')
    
    fax1.set_title('Distribution of diploid population')
    fax11.set_title('Distribution of tetraploid population')
    fax2.set_title('Distribution of resource')
    fax3.set_title('Population size of diploids')
    fax4.set_title('Population size of tetraploids')
    fax5.set_title('Amount of resources')

    faxcwD.set_title('Current weights of 2x')
    faxcwT.set_title('Current weights of 4x')
    faxcwDh.set_title('Current weights distribution of 2x')
    faxcwTh.set_title('Current weights distribution of 4x')
    faxerD.set_title('Current energy reserve of 2x')
    faxerT.set_title('Current energy reserve of 4x')
    faxerDh.set_title('2x energy reserve distribution')
    faxerTh.set_title('4x energy reserve distribution')

    faxProp.set_title('Proportion of tetraploids in population')

    faxmaturityDi.set_title('Average maturity time of 2x')
    faxmaturityTetra.set_title('Average maturity time of 4x')
    faxmaturityDih.set_title('2x maturity time distribution')
    faxmaturityTetrah.set_title('4x maturity time distribution')

    folder = 'results' #insert real folder name and parameters combination
    params = 'season100days_seeds0.5_10init10resby1_NND_per_0.099_alfa2_beta82_' 

    output_file =  open('%s/%s_%s.txt' % (folder, params, iter), 'w')
    output_file1 = open('%s/%s_indivs_%s.txt' % (folder, params, iter), 'w')
    output_file2 = open('%s/%s_parent_%s.txt' % (folder, params, iter), 'w')
    output_file1.write('%20s'%'Time' + '\t'+ \
                           '%20s'%'x' + '\t'+ \
                           '%20s'%'y' + '\t'+ \
                           '%20s'%'CurrentWeight'+'\t'+ \
                           '%20s'%'CurrentEnergyReserve'+'\t'+ \
                           '%20s'%'Ploidy'+'\t'+ \
                           '%20s'% 'Birthday'+'\t'+ \
                           '%20s'% 'Parent'+'\t'+ \
                           '%20s'% 'ResourcesInCell'+'\n')

    output_file3 = open('%s/%s_successfulrepr_%s.txt' % (folder, params, iter), 'w')
    output_file3.write('%20s'%'Time' + '\t'+ \
                           '%20s'%'x' + '\t'+ \
                           '%20s'%'y' + '\t'+ \
                           '%20s'%'CurrentWeight'+'\t'+ \
                           '%20s'%'CurrentEnergyReserve'+'\t'+ \
                           '%20s'%'Ploidy'+'\t'+ \
                           '%20s'% 'Birthday'+'\t'+ \
                           '%20s'% 'Parent'+'\t'+ \
                           '%20s'% 'ResourcesInCell'+'\n')

    output_file4 = open('%s/%s_main_indicators_%s.txt' % (folder, params, iter), 'a')
    output_file4.write('%20s'%'4x from 2x' + '\t'+ \
                           '%20s'%'4x from 4x' + '\t'+ \
                           '%20s'% '4x proportion'+'\n')
    output_file5 = open('%s/%s_maturity_weightreprDi_%s.txt' % (folder, params, iter), 'a')
    output_file5.write('%20s'%'timer' + '\t'+ \
                           '%20s'%'maturity time' + '\t'+ \
                           '%20s'% 'weight at reproduction'+'\n')
    output_file6 = open('%s/%s_maturity_weightreprTetra_%s.txt' % (folder, params, iter), 'a')
    output_file6.write('%20s'%'timer' + '\t'+ \
                           '%20s'%'maturity time' + '\t'+ \
                           '%20s'% 'weight at reproduction'+'\n')
    output_file7 = open('%s/%s_newborns_di_%s.txt' % (folder, params, iter), 'w')
    output_file7.write('%20s'%'Time' + '\t'+ \
                           '%20s'%'x' + '\t'+ \
                           '%20s'%'y' + '\t'+ \
                           '%20s'%'CurrentWeight'+'\t'+ \
                           '%20s'%'Ploidy'+'\n')
    output_file8 = open('%s/%s_newborns_tetra_%s.txt' % (folder, params, iter), 'w')
    output_file8.write('%20s'%'Time' + '\t'+ \
                           '%20s'%'x' + '\t'+ \
                           '%20s'%'y' + '\t'+ \
                           '%20s'%'CurrentWeight'+'\n')
    output_file9 = open('%s/%s_betadist_draws_%s.txt' % (folder, params, iter), 'w')
    output_file9.write('%20s'%'Time' + '\t'+ \
                           '%20s'%'x' + '\t'+ \
                           '%20s'%'y' + '\t'+ \
                           '%20s'%'Polyploid formation prob from Beta d' + '\t'+ \
                           '%20s'%'Polyploid formation prob from uniform d'+'\n')
    
    for timer in range(10000):
        meta.a_day_in_the_life(timer)
    output_file.close()
    output_file1.close()
    output_file2.close()
    output_file3.close()
    output_file4.close()
    output_file5.close()
    output_file6.close()
    output_file7.close()
    output_file8.close()
    output_file9.close()

    # saving figures 
    fig1.savefig('%s/%s_%s.png' % (folder, params, iter), dpi=300)
    fig2.savefig('%s/%s_weights_%s.png' % (folder, params, iter), dpi=300)
    fig4.savefig('%s/%s_energy_%s.png' % (folder, params, iter), dpi=300)
    fig3.savefig('%s/%s_PROPORTION_%s.png' % (folder, params, iter), dpi=300)
    fig5.savefig('%s/%s_parental_bars_%s.png' % (folder, params, iter), dpi=300)
    fig6.savefig('%s/%s_maturitytime_%s.png' % (folder, params, iter), dpi=300)
    fig7.savefig('%s/%s_maturitytimedist_%s.png' % (folder, params, iter), dpi=300)
    fig8.savefig('%s/%s_weightatrepr_%s.png' % (folder, params, iter), dpi=300)
    plt.close('all') 

    iter += 1  