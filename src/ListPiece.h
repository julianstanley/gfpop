#ifndef LISTPIECE_H
#define LISTPIECE_H

#include "Piece.h"
#include "Edge.h"
#include "ExternFunctions.h"

#include <math.h>

class ListPiece
{
private:
  Piece* head;
  Piece* currentPiece;
  Piece* lastPiece;

public:
  ListPiece();
  ~ListPiece();

  ///////  Simple list operations  ///////
  void setUniquePieceCostToInfinity();
  void setNewBounds(Interval newBounds);

  void reset();
  void copy(ListPiece  const& LP_edge);

  void reverse();

  void move();
  void initializeCurrentPiece();
  void initializeHeadWithFirstPoint(Point const& pt);

  void addCurrentPiecePlus1NotMove(Piece* newPiece);
  void addFirstPiece(Piece* newPiece);

  void shift(double parameter);
  void expDecay(double gamma);

  ///////  3 OPERATIONS in GFPOP ///////
  void LP_edges_constraint(ListPiece const& LP_state, Edge const& edge);
  void LP_edges_addPointAndPenaltyAndParentEdge(Point const& pt, Edge const& edge, unsigned int parentEdge);
  void LP_ts_Minimization(ListPiece& LP_edge, int parentEdge);

  ///////  operators up and down ///////
  void operatorUp(ListPiece const& LP_edge);
  void operatorDw(ListPiece const& LP_edge);

  void operatorSum(ListPiece& LP1, ListPiece& LP2);

  ///////  get info ///////
  void get_min_argmin_edge_ListPiece(double* response);
  void get_min_argmin_edge_ListPieceConstrained(double* response, Edge const& edge, bool forced);

  void show() const;
  void test();


};

#endif // LISTPIECE_H
