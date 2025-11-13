#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

#include "finite_element.hpp"

namespace fem {
template <typename T>
class Mesh {
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

  Mesh(Nodes const &nodes)
      : nodes(nodes),
        stiffnessMatrix(
            nodes.size() * FiniteElement::nNodeDofs,
            nodes.size() * FiniteElement::nNodeDofs),
        forceVector(nodes.size() * FiniteElement::nNodeDofs),
        displacementVector(nodes.size() * FiniteElement::nNodeDofs) {
    stiffnessMatrix.setZero();
    forceVector.setZero();
    displacementVector.setZero();
  }

  void printStiffnessMatrix() const {
    std::cout << stiffnessMatrix << std::endl;
  }

  void saveStiffnessMatrixToFile(const std::string &filename) const {
    std::ofstream file(filename);
    if (file.is_open()) {
      file << "Stiffness Matrix (" << stiffnessMatrix.rows() << "x"
           << stiffnessMatrix.cols() << "):\n";
      file << stiffnessMatrix << std::endl;
      file.close();
      std::cout << "Stiffness matrix saved to " << filename << std::endl;
    } else {
      std::cerr << "Error: Could not open " << filename << " for writing!"
                << std::endl;
    }
  }

  const Matrix &getStiffnessMatrix() const { return stiffnessMatrix; }

  // В mesh.hpp, метод calculateStiffnessMatrix:
  void calculateStiffnessMatrix(
      Value const &elasticityModulus, Value const &poissonRatio,
      typename FiniteElement::DifferentiationMatrix &dm) {
    typename FiniteElement::Nodes feNodes;
    for (Size i = 0; i < feNodes.size(); ++i) {
      feNodes(i) = nodes(i).coords;
    }

    typename FiniteElement::StiffnessMatrix sm;
    fe.integrateSubfields(sm, feNodes, elasticityModulus, poissonRatio, dm);

    stiffnessMatrix = sm;

    // ДОБАВЬТЕ ЭТО: пересчитайте dm для центра элемента
    // Найдём размеры элемента
    typename FiniteElement::Coordinates lsX, lsY, od;
    Value minx = feNodes(0)(0), maxx = feNodes(0)(0);
    Value miny = feNodes(0)(1), maxy = feNodes(0)(1);
    for (Size i = 0; i < feNodes.size(); ++i) {
      if (feNodes(i)(0) < minx) minx = feNodes(i)(0);
      if (feNodes(i)(0) > maxx) maxx = feNodes(i)(0);
      if (feNodes(i)(1) < miny) miny = feNodes(i)(1);
      if (feNodes(i)(1) > maxy) maxy = feNodes(i)(1);
    }
    Value a = maxx - minx;
    Value b = maxy - miny;

    // Вычисляем dm в центре элемента (xi=0, eta=0)
    fe.calculateDifferentiationMatrix(dm, 0.0, 0.0, a, b);
    this->dm = dm;
  }

  void calculateForceVector() {
    forceVector.setZero();

    for (Size i = 0; i < nodes.size(); ++i) {
      auto &forces = nodes(i).forces;
      if (forces.size() != 0) {
        for (Size j = 0; j < forces.size(); ++j) {
          Size dof = i * FiniteElement::nNodeDofs + forces(j).direction;
          forceVector(dof) = forces(j).value;
        }
      }
    }
  }

  void calculateDisplacementVector() {
    Matrix K = stiffnessMatrix;
    Vector F = forceVector;

    for (Size i = 0; i < nodes.size(); ++i) {
      auto &disps = nodes(i).disps;
      if (disps.size() != 0) {
        for (Size j = 0; j < disps.size(); ++j) {
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

  void writeParaViewVtk() {
    std::ofstream vtk("output.vtk");
    if (!vtk.is_open()) {
      std::cerr << "Error: Could not open output.vtk for writing!" << std::endl;
      return;
    }
    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "Finite Element Solution\n";
    vtk << "ASCII\n";
    vtk << "DATASET UNSTRUCTURED_GRID\n\n";

    vtk << "POINTS " << nodes.size() << " float\n";
    for (Size i = 0; i < nodes.size(); ++i) {
      vtk << nodes(i).coords(0) << " " << nodes(i).coords(1) << " 0.0\n";
    }

    vtk << "\nCELLS " << nElems << " " << nElems * 5 << "\n";
    vtk << "4 0 1 2 3\n";

    vtk << "\nCELL_TYPES " << nElems << "\n";
    vtk << "9\n";

    vtk << "\nPOINT_DATA " << nodes.size() << "\n";
    vtk << "VECTORS displacement float\n";
    for (Size i = 0; i < nodes.size(); ++i) {
      vtk << displacementVector(FiniteElement::nNodeDofs * i) << " "
          << displacementVector(FiniteElement::nNodeDofs * i + 1) << " 0.0\n";
    }
    vtk.close();
    std::cout << "Successfully wrote output.vtk" << std::endl;
  }

  const Vector &getDisplacementVector() const { return displacementVector; }

  const Vector getNDS(Value const &elastMod) {
    return fe.calculateNDS(displacementVector, elastMod, this->dm);
  }

 private:
  Nodes nodes;
  FiniteElement fe;
  typename FiniteElement::DifferentiationMatrix dm;
  Matrix stiffnessMatrix;
  Vector forceVector, displacementVector;
};
}  // namespace fem