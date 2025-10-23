#include "mesh.hpp"
#include <iostream>
#include <fstream>

int main()
{
    using Value = float;
    using Mesh = fem::Mesh<Value>;

    // typename Mesh::Nodes nodes {{{0, 0}, {}, {{0, 0}, {1, 0}}}, // x - 0 ||| y - 1
    //                             {{1, 0}, {}, {{0, 0}, {1, 0}}}, 
    //                             {{1, 1}, {{1, -1000}}, {}}, 
    //                             {{0, 1}, {{1, -1000}}, {}}};
    typename Mesh::Nodes nodes(4);
    nodes(0).coords = {0, 0};
    nodes(1).coords = {1, 0};
    nodes(2).coords = {1, 1};
    nodes(3).coords = {0, 1};

    nodes(0).disps.resize(2);
    nodes(0).disps(0) = {0, 0}; 
    nodes(0).disps(1) = {1, 0};     
    nodes(1).disps.resize(2);
    nodes(1).disps(0) = {0, 0}; 
    nodes(1).disps(1) = {1, 0}; 

    nodes(2).forces.resize(1);
    nodes(2).forces(0) = {1, -1000}; 
    nodes(3).forces.resize(1);
    nodes(3).forces(0) = {1, -1000}; 



            
    Mesh mesh(nodes);
    Value const elastMod = 200000;
    Value const poissRat = 0.3;

    std::cout << "Starting FEM calculation..." << std::endl;

    mesh.calculateStiffnessMatrix(elastMod, poissRat);
    std::cout << "Stiffness matrix calculated." << std::endl;

    // ДОБАВЬТЕ ЭТО - вывод матрицы жесткости в консоль
    std::cout << "Stiffness matrix:" << std::endl;
    mesh.printStiffnessMatrix(); // Этот метод нужно добавить в класс Mesh

    // ДОБАВЬТЕ ЭТО - сохранение матрицы в файл
    mesh.saveStiffnessMatrixToFile("stiffness_matrix.txt");

    mesh.calculateForceVector();
    std::cout << "Force vector calculated." << std::endl;

    mesh.calculateDisplacementVector();
    std::cout << "Displacement vector calculated." << std::endl;

    mesh.writeParaViewVtk();
    std::cout << "VTK file written." << std::endl;

    std::cout << "Calculation completed. Check output.vtk and stiffness_matrix.txt" << std::endl;
    std::cout << "Press Enter to continue...";
    std::cin.get();

    return 0;
}