// to add Finite Element to ff++
class EConstantTypeOfFE :public E_F0
{ 
//  using namespace   Fem2D;
  Fem2D::TypeOfFE * v;
public:
  AnyType operator()(Stack ) const { /*cout << " ()" << v << endl*/;return SetAny<Fem2D::TypeOfFE*>(v);}
  EConstantTypeOfFE( Fem2D::TypeOfFE * o):v(o) { /*cout << "New constant " << o << endl;*/}
  size_t nbitem() const { assert(v);return v->N ;} 
   operator aType () const { return atype<Fem2D::TypeOfFE*>();} 
};


struct AddNewFE {  
  AddNewFE (const char * FEname,Fem2D::TypeOfFE* tfe) 
  {
    ffassert(tfe); // check 
    Global.New(FEname, Type_Expr(atype<Fem2D::TypeOfFE*>() ,new  EConstantTypeOfFE(tfe)));
  }
};