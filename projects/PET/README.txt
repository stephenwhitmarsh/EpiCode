Trials d�finis par : 
	- d�but 2 secondes apr�s la pointe d'une IED et fin 2 secondes avant la pointe de l'IED suivante
	- doit durer au moins 5 secondes
	- ne doit pas contenir de marqueur BAD
	- si trial sup�rieur � 10 seconde, re-coup� en plusieurs trials de 5 secondes


PREPARER
1/ V�rifier les marqueurs Muse avec PET_verify_markers.m
2/ Pr�parer le patient sur PET_setparams.m puis aller dans le dossier \\l2export\iss02.charpier\analyses\vn_pet\data\wrong_markers

SPIKE SORTING
3/ Ex�cuter PET_prepare_spykingcircus.m en indiquant le num�ro du patient

4/ Lancer spyking circus sur la machine virtuelle (ou X2Go si �a prend trop de m�moire sur l'ordi local)
4bis/ Ou, lancer Spyking-circus sur le cluser : 
	- ouvrir Git Bash
	- ssh login02 (connexion au cluster)
	- alloc_spykingcircus
	- ssh ... (� remplacer par le nom du node : hmb001, lmb043, etc.)
	- petdata (ouvrir les donn�es pour spyking circus)
	- cd ... (� remplacer par le nom du patient)
	- cd p1
	- cd ... (� remplacer par le nom du groupe d'�lectrode : mAmT2 par exemple)
	- module load spyking-circus/1.08
	- spyking-circus SpykingCircus.params -c 28

5/ Trier sur Phy entre noise/mu/god, et merge certains clusters si besoin


ANALYSE MATLAB
6/ Ex�cuter PET_project_per_patient sur le cluster.
	- ouvrir Git Bash
	- ssh login02 (connexion au cluster)
	- pet (ouvrir le dossier avec les scripts)
	- sbatch --array=PAT_NR pet_project_per_patient_slurm.sh (PAT_NR : a remplacer par le numero du ou des patients (sans espace apr�s le "="): 1, 1-4, 1,3,7 etc.)
	- squeueme : voir si c'est encore en train de tourner ou bien si c'est termin�
	- aller dans le dossier \\l2export\iss02.charpier\analyses\vn_pet\slurm_error : voir s'il y a eu une erreur pendant l'ex�cution

7/ - Verifier les plots et si besoin revenir sur Phy pour modifier la classification des units
   - Puis relancer les scripts d'analyse

8/ Une fois toutes les analyse faites, avec le plus de patients possible, lancer : PET_classification_celltype. 
Ne pas utiliser les r�sultats si les plots ne sont pas convaincants (ne montrent pas une distribution bimodale)