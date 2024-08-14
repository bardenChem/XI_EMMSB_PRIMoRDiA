#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#=======================================================================
'''
Autor: Dr. Igor Barden Grillo
email: barden.igor@gmail.com
Minicurso: Métodos químico-quânticos e descritores de 
reatividade aplicados a sistemas biológicos

Script desenvolvido para rodar as simulações de caminho de reação enzimá-
ca durante o mini-curso da XI EMMSB do LNCC.  
'''

#=======================================================================

from pDynamo_Scripts import Scripts # importa a biblioteca principal do pDynamo_Scripts
import SimulationSystem             #
import os, sys                      # bibliotecas auxiliares comumns do Python     

#==============================================
def Prepare_MM_System():
	'''
	Função para configurar o sistema a partir de arquivos de campo de 
	força do AMBER.  
	'''
	#Estrutura de dicionário com as opções e os respectivos parâmetros
	_opt_pars = {
		#Tipo de input : arquivo de campo de força do Amber 
		"Input_Type":"amber"                     , 
		#Tipo de simulação: otimização de geomteria             
		"simulation_type":"Geometry_Optimization", 
		#Formato da trajetória para ser salva
		"save_format":".dcd",
		#Frequência com que as estruturas são salvas
		"save_frequency":20,
		#Tolerância de variação da estrutura para parada da simulação
		"rmsGradient":2.0,
		#Número máximo de iterações
		"maxIterations":2200,
		#Nome do arquivo de trajetória para ser salvo
		"trajectory_name":"opt_full_tim.ptGeo",
		#Nome do arquivo com as coordenadas
		"crd_file":"7tim.crd",
		#Nome do arquivo com as coordenadas
		"top_file":"7tim.top"
	}

	#Objeto da classe "scripts" com o nome da pasta para ser gerada
	test_03 = Scripts("Prep_system")
	#Método de classe para configurar o sistema
	test_03.Set_System(_opt_pars)
	#Método de classe para rodar a simulação 
	test_03.Run_Simulation(_opt_pars)
	#Método de classe para salvar o sistema em PDB e em PKL (arquivo serializado do pDynamo)
	test_03.SaveSystem()

#-----------------------------------------------
def Prepare_Prune_System():
	'''
	'''
	#Estrutura de dicionário com as opções e os respectivos parâmetros
	_prune_pars= {
		#Tipo de input : arquivo serializado do pDynamo
		"Input_Type":"pkl",
		#Nome do arquivo PKL
		"pkl_file":"Prep_system/7tim.pkl",
		#Tipo de simulação
		"simulation_type":"Geometry_Optimization",
		#Sinalização para fazer o recorte esférico a partir do átomo C02 do ligante/substrato
		"spherical_prune":"*:LIG.248:C02",
		#Raio de corte do sistema 
		"spherical_prune_radius":25.0,
		#Sinalizção para congelar o movimento dos átomos a partir de certo raio em relação ao átomo C02 do ligante/substrato
		"set_fixed_atoms":"*:LIG.248:C02",
		#Raio para definir a partir de que distância os átomos serão fixados
		"free_atoms_radius":20.0,
		#Formato para salvar a trajetória de otimização de geometria
		"save_format":".dcd",
		#Frequência com que as estruturas são salvas
		"save_frequency":20,
		#Tolerância de variação da estrutura para parada da simulação
		"rmsGradient":1,
		#Número máximo de iterações
		"maxIterations":2200,
		#Nome do arquivo de trajetória para ser salvo
		"trajectory_name":"opt_pruned_tim.ptGeo"
	}
	#Objeto da classe "scripts" com o nome da pasta para ser gerada
	test_03_b = Scripts("Prep_prune")
	#Método de classe para configurar o sistema
	test_03_b.Set_System(_prune_pars)
	#Método de classe para rodar a simulação 
	test_03_b.Run_Simulation(_prune_pars)
	#Método de classe para salvar o sistema em PDB e em PKL (arquivo serializado do pDynamo)
	test_03_b.SaveSystem("7tim_optMM")


#-----------------------------------------------
def Set_QC_MM(_hamiltonian):
	'''
	'''
	#Estrutura de dicionário com as opções e os respectivos parâmetros
	_qc_mmpars = {
		"Input_Type":"pkl",
		"pkl_file":"Prep_prune/7tim_optMM.pkl",
		"set_energy_model":"QM",
		"Hamiltonian":_hamiltonian,
		"method_class":"SMO",
		"set_qc_region":"yes",
		"residue_patterns":["*:LIG.248:*","*:GLU.164:*","*:HIE.94:*"],
		"QCcharge":-3,
		"save_format":".dcd",
		"save_frequency":20,
		"rmsGradient":0.1,
		"maxIterations":2200,
		"trajectory_name":"opt_qcmm_tim_"+_hamiltonian+".ptGeo",
		"simulation_type":"Geometry_Optimization"
	}

	test_03_c = Scripts("QM_opt")
	test_03_c.Set_System(_qc_mmpars)
	test_03_c.Run_Simulation(_qc_mmpars)
	test_03_c.SaveSystem("7tim_"+_hamiltonian+"_opt_PF")


#-------------------------------------------------
def Multiple_Distance(_hamiltonian):
	'''
	'''
	#Estrutura de dicionário com as opções e os respectivos parâmetros
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("QM_opt","7tim_"+_hamiltonian+"_opt_PF.pkl"),		
		"set_reaction_crd":1,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		#"atoms_rc2":["*:LIG.*:O06","*:HIE.94:HE2","*:HIE.94:NE2"],
		"type":"multipleDistance",
		"maxIterations":2200,
		"mass_constraint":"True",
	}
	#Estrutura de dicionário com as opções e os respectivos parâmetros
	scan1_parameters = {
		"simulation_type":"Relaxed_Surface_Scan",
		"dincre_rc1":0.1,
		"optmizer":"ConjugatedGradient",
		"maxIterations":2200,
		"nsteps_rc1":20,
		"force_constants":[4000.0,4000.0]

	}
	#test simple distance
	test_01 = Scripts("Multiple_Distance_"+_hamiltonian)
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(scan1_parameters)
	test_01.SaveSystem("Multiple_DistanceScan")


#=======================================================================
def Run_Mopac_Refinement():
	'''
	'''
	#Estrutura de dicionário com as opções e os respectivos parâmetros
	system_parameters = {
		"Input_Type":"pkl",		
		"pkl_file":os.path.join("QM_opt","7tim_rm1_opt_PF.pkl"),
		"set_qc_region":"yes",
		"set_energy_model":"QM",
		"Hamiltonian":"rm1",
		"method_class":"SMO",
		"center_atom":"*:LIG.*:H02",
		"radius": 5.0,
		"set_reaction_crd":1,
		"atoms_rc1":["*:LIG.*:C02","*:LIG.*:H02","*:GLU.164:OE2"],
		"type":"Distance",
		"mass_constraint":"True",
	}

	_path   = "Multiple_Distance_rm1/ScanTraj.ptGeo"
	methods = ["am1","pm3","rm1","pm6","pm7"]
	#Estrutura de dicionário com as opções e os respectivos parâmetros
	simulation_parameters = { "xnbins":20			    ,
				   "source_folder":_path                , 
				   "folder":"Mopac_Refinement"          ,
				   "charge":-2		                    ,
				   "multiplicity":1 	                ,
				   "methods_lists":methods              ,   
				   "NmaxThreads": 4                     ,
				   "simulation_type":"Energy_Refinement",
				   "Software":"mopac"	}				  
				
	#------------------------------------
	test_01 = Scripts("Mopac_Refinement")
	test_01.Set_System(system_parameters)
	test_01.Run_Simulation(simulation_parameters)
	test_01.SaveSystem()
	#-----------------------------------
#-----------------------------------------------
def Run_Test():
	'''
	Test 
	'''
	
	#Confere de o caminho "Prep_system/7tim.pkl",
	# se não existir executa a preparação do sistema
	if not os.path.exists( "Prep_system/7tim.pkl" ):
		print("Preparing System\n")
		print("======================================")
		Prepare_MM_System()
	
	#Confere de o caminho "Prep_prune/7tim_optMM.pkl",
	# se não existir executa a preparação do sistema
	if not os.path.exists("Prep_prune/7tim_optMM.pkl"):
		print("Pruning System\n")
		print("======================================")
		Prepare_Prune_System()
	
	#Confere de o caminho "QM_opt/7tim_rm1_opt_PF.pkl",
	# se não existir executa a preparação do sistema
	if not os.path.exists( "QM_opt/7tim_rm1_opt_PF.pkl"):
		print("QCMM Optimization")
		print("======================================")
		Set_QC_MM("rm1")

	#Executa a varredura relaxada da coordenada de reação.
	Multiple_Distance("rm1")
	#Executa o refinamento de energia no mopac 
	Run_Mopac_Refinement()
	
#===================================
if __name__ == '__main__': Run_Test()
