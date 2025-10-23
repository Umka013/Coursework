#pragma once

#include "finite_element.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

namespace fem
{
    template <typename T>
    class Mesh
    {
    public:
        using Value = T;
        using Size = unsigned long long int;
        using FiniteElement = fem::FiniteElement<Value>;
        Size const static nElems = 1;
        Size const static nMeshDofs = nElems * FiniteElement::nElemDofs;
        using Vector = Eigen::VectorX<Value>;
        using Matrix = Eigen::MatrixX<Value>;

        struct Cond {
            Size direction;
            Value value;
        };

        using Conds = Eigen::VectorX<Cond>;

        struct Node {
            typename FiniteElement::Coordinates coords;
            Conds forces, disps;
        };

        using Nodes = Eigen::VectorX<Node>;

        Mesh(Nodes const& nodes) : 
            nodes(nodes),
            stiffnessMatrix(nodes.size() * FiniteElement::nNodeDofs, nodes.size() * FiniteElement::nNodeDofs),
            forceVector(nodes.size() * FiniteElement::nNodeDofs),
            displacementVector(nodes.size() * FiniteElement::nNodeDofs)
        {
            stiffnessMatrix.setZero();
            forceVector.setZero();
            displacementVector.setZero();
        }

        void printStiffnessMatrix() const
        {
            std::cout << stiffnessMatrix << std::endl;
        }

        void saveStiffnessMatrixToFile(const std::string &filename) const
        {
            std::ofstream file(filename);
            if (file.is_open())
            {
                file << "Stiffness Matrix (" << stiffnessMatrix.rows()
                     << "x" << stiffnessMatrix.cols() << "):\n";
                file << stiffnessMatrix << std::endl;
                file.close();
                std::cout << "Stiffness matrix saved to " << filename << std::endl;
            }
            else
            {
                std::cerr << "Error: Could not open " << filename << " for writing!" << std::endl;
            }
        }

        const Matrix &getStiffnessMatrix() const
        {
            return stiffnessMatrix;
        }

        void calculateStiffnessMatrix(
            Value const &elasticityModulus,
            Value const &poissonRatio)
        {
            typename FiniteElement::Nodes feNodes;
            for(Size i = 0; i < feNodes.size(); ++i) {
                feNodes(i) = nodes(i).coords;
            }

            typename FiniteElement::StiffnessMatrix sm;
            fe.calculateStiffnessMatrix(sm, feNodes, elasticityModulus, poissonRatio);
            
            stiffnessMatrix = sm;
        }

        void calculateForceVector()
        {
            forceVector.setZero();
            
            for(Size i = 0; i < nodes.size(); ++i) {
                auto &forces = nodes(i).forces;
                if (forces.size() != 0) {
                    for(Size j = 0; j < forces.size(); ++j) {
                        Size dof = i * FiniteElement::nNodeDofs + forces(j).direction;
                        forceVector(dof) = forces(j).value;
                    }
                }
            }            
        }

        void calculateDisplacementVector()
        {
            Matrix K = stiffnessMatrix;
            Vector F = forceVector;

            for (Size i = 0; i < nodes.size(); ++i)
            {
                auto &disps = nodes(i).disps;
                if (disps.size() != 0) {
                    for(Size j = 0; j < disps.size(); ++j) {
                        Size dof = i * FiniteElement::nNodeDofs + disps(j).direction;
                        K.row(dof).setZero();
                        K.col(dof).setZero();
                        K(dof, dof) = 1.0;
                        F(dof) = disps(j).value;
                    }
                }
            }

            displacementVector = K.lu().solve(F);
        }

        void writeParaViewVtk()
        {
            std::ofstream vtk("output.vtk");
            if (!vtk.is_open())
            {
                std::cerr << "Error: Could not open output.vtk for writing!" << std::endl;
                return;
            }
            vtk << "# vtk DataFile Version 3.0\n";
            vtk << "Finite Element Solution\n";
            vtk << "ASCII\n";
            vtk << "DATASET UNSTRUCTURED_GRID\n\n";

            vtk << "POINTS " << nodes.size() << " float\n";
            for (Size i = 0; i < nodes.size(); ++i)
            {
                vtk << nodes(i).coords(0) << " " << nodes(i).coords(1) << " 0.0\n";
            }

            vtk << "\nCELLS " << nElems << " " << nElems * 5 << "\n";
            vtk << "4 0 1 2 3\n";

            vtk << "\nCELL_TYPES " << nElems << "\n";
            vtk << "9\n";

            vtk << "\nPOINT_DATA " << nodes.size() << "\n";
            vtk << "VECTORS displacement float\n";
            for (Size i = 0; i < nodes.size(); ++i)
            {
                vtk << displacementVector(FiniteElement::nNodeDofs * i) << " "
                    << displacementVector(FiniteElement::nNodeDofs * i + 1) << " 0.0\n";
            }
            vtk.close();
            std::cout << "Successfully wrote output.vtk" << std::endl;
        }

    private:
        Nodes nodes;
        FiniteElement fe;
        Matrix stiffnessMatrix;
        Vector forceVector, displacementVector;
    };
}