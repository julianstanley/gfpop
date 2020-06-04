#ifndef PIECE_H
#define PIECE_H

#include "Edge.h"
#include "Data.h"

#include "Cost.h"
#include "Interval.h"

#include<vector>
#include<string>

#include <fstream> ///write in a file

class Piece
{
  public:
    Cost m_cost;
    Interval m_interval;
    int m_parentEdge;
    Piece *nxt;   /// pointer to next piece

    Piece();
    Piece(Cost const& cost, Interval const& inter = Interval(), int myEdge = -1);
    Piece(const Piece* piece); ///COPY CONSTRUCTOR => copy only the first Piece. piece -> nxt = NULL
    ~Piece();
    Piece* copy();

    void addCostAndPenaltyAndParentEdge(Cost const& cost, double penalty, unsigned int parentEdge);

    ///
    ///
    Interval intervalMinLessUp(double bound, double currentValue, bool constPiece);
    Interval intervalMinLessDw(double bound, double currentValue, bool constPiece);
    Piece* pastePieceUp(const Piece* NXTPiece, Interval const& decrInter);
    Piece* pastePieceDw(const Piece* NXTPiece, Interval const& decrInter);

    Piece* pieceGenerator(Piece* Q1, Piece* Q2, int Bound_Q2_Minus_Q1, double M);
    Piece* piece0(Piece* Q1, Piece* Q2, Interval interToPaste, int& Q2_Minus_Q1);
    Piece* piece1(Piece* Q1, Piece* Q2, Interval interToPaste, Interval interRoots, int& Q2_Minus_Q1);
    Piece* piece2(Piece* Q1, Piece* Q2, Interval interToPaste, Interval interRoots, int& Q2_Minus_Q1);

    void get_min_argmin_edge(double* response);
    ///
    ///

    void show();

};

std::ostream &operator>>(std::ostream &flux, Piece* piece);


#endif // PIECE_H
