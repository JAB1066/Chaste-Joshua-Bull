# CMake generated Testfile for 
# Source directory: /scratch/eclipse/workspace/chaste/projects/JoshuaBull/test
# Build directory: /scratch/eclipse/workspace/chaste-build/projects/JoshuaBull/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TestHello "/scratch/eclipse/workspace/chaste-build/projects/JoshuaBull/test/TestHelloRunner")
set_tests_properties(TestHello PROPERTIES  LABELS "project_JoshuaBull;Continuous" PROCESSORS "1" WORKING_DIRECTORY "/scratch/eclipse/workspace/chaste/")
add_test(TestRunningTumourSpheroidSimulationsTutorialEdited "/scratch/eclipse/workspace/chaste-build/projects/JoshuaBull/test/TestRunningTumourSpheroidSimulationsTutorialEditedRunner")
set_tests_properties(TestRunningTumourSpheroidSimulationsTutorialEdited PROPERTIES  LABELS "project_JoshuaBull;Continuous" PROCESSORS "1" WORKING_DIRECTORY "/scratch/eclipse/workspace/chaste/")
