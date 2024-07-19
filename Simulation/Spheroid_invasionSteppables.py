from cc3d.core.PySteppables import *

import numpy as np
import os
import csv 

front_angle = 3.14/4
F_ECM = 1.0
thresh = {{thresh}} # 0.05 or 5% threshold for the hill function
nThresh = {{nThresh}} 

f_R = {{f_R}}
PL0 = {{PL0}}#This one is used with the active ECM feedback mechanism
PF = 1.0

cLF = {{cLF}}

xDim = {{xDim}}

totalSimTime = {{totalSimTime}}
measuresDataFreq = {{measuresDataFreq}}
cellDataFreq = {{cellDataFreq}}

class CalculateARandFrontLeaders(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)
        
        
    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        x_min = 0.0
        x_max = 0.0
        y_min = 0.0
        y_max = 0.0
        z_min = 0.0
        z_max = 0.0
        
        org_com = np.array([0.0, 0.0, 0.0])
        fol_com = np.array([0.0, 0.0, 0.0])
        leader_com = np.array([0.0, 0.0, 0.0])
        
        front_count = 0
        cellPos = np.zeros((len(self.cell_list), 3))
        for num, cell in enumerate(self.cell_list):
            cellPos[num, :] = np.array([cell.xCOM, cell.yCOM, cell.zCOM])
            
            org_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
            
            if cell.type == self.cell_type.FOLLOWER:
                fol_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
                
            if cell.type == self.cell_type.LEADER:
                leader_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])            
                
        org_com = org_com/len(self.cell_list)
        fol_com = fol_com/len(self.cell_list_by_type(self.cell_type.FOLLOWER))
        leader_com = leader_com/len(self.cell_list_by_type(self.cell_type.LEADER))
        
        self.shared_steppable_vars['current_org_com'] = org_com
        self.shared_steppable_vars['current_fol_com'] = fol_com
        self.shared_steppable_vars['current_leader_com'] = leader_com
        
        max_values = np.amax(cellPos, axis=0)
        min_values = np.amin(cellPos, axis=0)
        xMax, yMax, zMax = max_values
        xMin, yMin, zMin = min_values

        if (yMax-yMin) > (xMax-xMin): # Check Y axis is larger than x-axis
            self.shared_steppable_vars['org_AR'] = (xMax-xMin)/(yMax-yMin)
        else: # X-axis is larger than y-axis
            self.shared_steppable_vars['org_AR'] = (yMax-yMin)/(xMax-xMin)

        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            vec = -fol_com + np.array([cell.xCOM, cell.yCOM, cell.zCOM])
            x_unit_vector = np.array([1,0,0])
            dot_prod = np.dot(vec,x_unit_vector)
            cos_theta = dot_prod/self.vector_norm(vec)
            
            if (dot_prod > 40*np.cos(front_angle) ):#40*cos(theta) defines our conical region
                front_count += 1

        self.shared_steppable_vars['front_leaders'] = front_count

        if mcs%50 == 0:
            print('FINISHED STEP ' + str(mcs))

class Calculate_P(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)
        
    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        initial_org_com = np.array([0.0, 0.0, 0.0])
        initial_L_com = np.array([0.0, 0.0, 0.0])
        initial_fol_com = np.array([0.0, 0.0, 0.0])
        self.shared_steppable_vars['initial_org_com'] = np.array([0.0, 0.0, 0.0])
        self.shared_steppable_vars['initial_leader_com'] = np.array([0.0, 0.0, 0.0])
        self.shared_steppable_vars['initial_fol_com'] = np.array([0.0, 0.0, 0.0])

        for cell in self.cell_list:
            initial_org_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        initial_org_com /= len(self.cell_list)
        self.shared_steppable_vars['initial_org_com'] = initial_org_com
        
        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            vec = np.array([cell.xCOM,cell.yCOM,cell.zCOM]) - initial_org_com
            if self.vector_norm(vec) > 23:
                cell.type = self.cell_type.FOLLOWER

        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            initial_L_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        initial_L_com /= len(self.cell_list_by_type(self.cell_type.LEADER))
        self.shared_steppable_vars['initial_leader_com'] = initial_L_com

        for cell in self.cell_list_by_type(self.cell_type.FOLLOWER):
            initial_fol_com += np.array([cell.xCOM, cell.yCOM, cell.zCOM])
        initial_fol_com /= len(self.cell_list_by_type(self.cell_type.FOLLOWER))
        self.shared_steppable_vars['initial_fol_com'] = initial_fol_com
                        

        print ("percentage of leaders", 100*len(self.cell_list_by_type(self.cell_type.LEADER) )/len(self.cell_list) )
        print ("num of leaders", len(self.cell_list_by_type(self.cell_type.LEADER) )  )

        for cell in self.cell_list:

            
            if cell.type == self.cell_type.LEADER:
                # cell.dict['PL0'] = 10.0
                cell.dict['beta'] = 0.0
                
            if cell.type == self.cell_type.FOLLOWER:
                cell.dict['P'] = 0.0
##########################################################################################################################    

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """
        # CCCP = 1.0 Cell-Cell Communication Parameter
        # delta_t = 1/60*90/1 1 MCS = 1 minute
        global PL0
        global f_R
        tot_leaders = len(self.cell_list_by_type(self.cell_type.LEADER) )

        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            free_boundary_norm = 0.0
            flag_boundary = 0
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if not neighbor:
                    fbn = common_surface_area/cell.surface #fbn =  free boundary normalized
                    cell.dict['beta'] = fbn**nThresh/(fbn**nThresh + thresh**nThresh)
                    flag_boundary = 1
            if flag_boundary == 0:
                cell.dict['beta'] = 0.0

        n_activated_leader = 0

        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            if cell.dict['beta'] > (1.0-thresh):
                n_activated_leader += 1
        fr_net = 0.0
        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            beta = cell.dict['beta']
            fr_net = n_activated_leader/len(self.cell_list)*beta*f_R
            if beta == 0.0:
                P = PL0 
            else:
                # P = PL0*(1.0 - 2.71**(-fr_net/F_ECM) )
                P = PL0*(1.0 - np.exp(-fr_net/F_ECM) )

            cell.dict['P'] = P


        PF = 0.0
        beta_tot = 0.0
        for cell in self.cell_list_by_type(self.cell_type.LEADER):
            beta = cell.dict['beta']
            PF += cell.dict['P']*beta
            beta_tot += beta
        if beta_tot != 0.0:
            PF = PF/beta_tot

        for cell in self.cell_list_by_type(self.cell_type.FOLLOWER):
            cell.dict['P'] = cLF * PF # vary cLF from 0.1 to 1

        currentDir = self.output_dir
        splitPath = os.path.split(currentDir)
        fileName = 'measures_f' + str(measuresDataFreq) +'.csv'

        output_path = os.path.join(splitPath[0], fileName)

        if  not os.path.isfile(output_path):  # false--> true
            # print('CellInfo.csv DOES NOT EXIST. IT WILL BE CREATED.')
            with open(output_path, 'a') as fout:
                writer = csv.writer(fout, delimiter=',')
                writer.writerow(['mcs', 'orgAvgSpeed', 'leaderAvgSpeed', 'followerAvgSpeed', 'orgAR', 'frontLeadPerc'])
                writer.writerow([mcs, 0, 0, 0, self.shared_steppable_vars['org_AR'], self.shared_steppable_vars['front_leaders']/tot_leaders])

        else:
            if mcs%measuresDataFreq == 0:
                delta_org = self.shared_steppable_vars['current_org_com'] - self.shared_steppable_vars['initial_org_com'] # NOT INSTANTANEOUS, BUT AVG
                delta_leader = self.shared_steppable_vars['current_leader_com'] - self.shared_steppable_vars['initial_leader_com'] # NOT INSTANTANEOUS, BUT AVG
                delta_fol = self.shared_steppable_vars['current_fol_com'] - self.shared_steppable_vars['initial_fol_com'] # NOT INSTANTANEOUS, BUT AVG
            
                with open(output_path, 'a') as fout:
                    writer = csv.writer(fout, delimiter=',')
                    writer.writerow([mcs, (delta_org[0])/mcs, (delta_leader[0])/mcs, (delta_fol[0])/mcs, self.shared_steppable_vars['org_AR'], self.shared_steppable_vars['front_leaders']/tot_leaders])


    def finish(self):
        """
        Finish Function is called after the last MCS
        """
        pass


class Spheroid_invasionSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
#         initial_pos_population = np.array([0.0,0.0,0.0])
        for cell in self.cell_list:
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 1.0
            cell.targetSurface = 1.1*cell.surface
            cell.lambdaSurface = 0.1

    def step(self,mcs):
        """
        type here the code that will run every frequency MCS
        :param mcs: current Monte Carlo step
        """

        for cell in self.cell_list:
            cell.dict["migDir"] = np.array([[1, 0, 0]]).T

        for cell in self.cell_list_by_type(self.cell_type.LEADER):

            avgPL += cell.dict['P']

            # For ECM feedback
            muMag = cell.dict['P']
            cell.lambdaVecX = -1.0*muMag*cell.dict['migDir'][0,0]
            cell.lambdaVecY = -1.0*muMag*cell.dict['migDir'][1,0]
            cell.lambdaVecZ = -1.0*muMag*cell.dict['migDir'][2,0]

        avgPL /= len(self.cell_list_by_type(self.cell_type.LEADER))

        for cell in self.cell_list_by_type(self.cell_type.FOLLOWER):
            muMag = ( (cell.dict['P'] <=1)*(1 - cell.dict['P']) + cell.dict['P'] )
            
            cell.lambdaVecX = -1.0*muMag*cell.dict['migDir'][0,0]
            cell.lambdaVecY = -1.0*muMag*cell.dict['migDir'][1,0]
            cell.lambdaVecZ = -1.0*muMag*cell.dict['migDir'][2,0]

            # First attempt at changes to follower-leader comms
            ## cell.lambdaVecX = -1 * ((cell.dict['P'] <=1) * (cell.dict['P'] + cLF*(avgPL - cell.dict['P'])) )

            # For no ECM feedback
            # cell.lambdaVecX = -1.0
            

        currentDir = self.output_dir
        splitPath = os.path.split(currentDir)
        fileName1 = 'cellInfo' + '.csv'

        output_path1 = os.path.join(splitPath[0], fileName1)

        if  not os.path.isfile(output_path1):  # false--> true
            # print('CellInfo.csv DOES NOT EXIST. IT WILL BE CREATED.')
            with open(output_path1, 'a') as fout:
                writer = csv.writer(fout, delimiter=',')
                writer.writerow(['cellID', 'cell.xCOM', 'cell.yCOM', 'cell.zCOM', 'axes0', 'axes1', 'axes2', 'surface', 'volume','cellType', 'cell.LambdaVecX', 'cell.LambdaVecY', 'cell.LambdaVecZ', 'cell.P'])
                for cell in self.cell_list:
                    axes=self.momentOfInertiaPlugin.getSemiaxes(cell)
                    writer.writerow([cell.id, cell.xCOM, cell.yCOM, cell.zCOM, axes[0], axes[1], axes[2], cell.surface, cell.volume, cell.type, cell.lambdaVecX, cell.lambdaVecY, cell.lambdaVecZ, cell.dict['P']])
        else:
            if mcs%cellDataFreq == 0:
                with open(output_path1, 'a') as fout:
                    writer = csv.writer(fout, delimiter=',')
                    for cell in self.cell_list:
                        axes=self.momentOfInertiaPlugin.getSemiaxes(cell)
                        writer.writerow([cell.id, cell.xCOM, cell.yCOM, cell.zCOM, axes[0], axes[1], axes[2], cell.surface, cell.volume, cell.type, cell.lambdaVecX, cell.lambdaVecY, cell.lambdaVecZ, cell.dict['P']])


        if mcs >= totalSimTime:
            self.stop_simulation()

    def finish(self):
        """
        Finish Function is called after the last MCS
        """