#ifndef INCLUDED_actor_PMActor_HH
#define INCLUDED_actor_PMActor_HH

#include <Eigen/Dense>
#include <scheme/io/dump_pdb_atom.hh>

#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzTransform.hh>

namespace scheme {
namespace actor {



	template<class _Position>
	struct PMActor {
		//typedef ::Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;
		///@brief Position type, leave out to make actor "Fixed"
		typedef _Position Position;
		typedef typename Position::Scalar Float;
		// typedef EigenXform Position;
		// typedef typename Position::Scalar Float;
		typedef PMActor THIS;
		typedef Eigen::Matrix<Float,3,1> V3;

		Position position_;
		char aa_,ss_;
		int index_;

		PMActor() : position_(Position::Identity()), aa_('-'), ss_('-'), index_(0) {}

		PMActor(Position const & p, char aa='-', char ss='-', int i=0) : position_(p), aa_(aa), ss_(ss),index_(i) {}

		// SCPMActor( 
		// 	core::conformation::Residue const & res,
		// 	char aa='-', 
		// 	char ss='-', 
		// 	int i=0
		// )  : aa_(aa), ss_(ss), index_(i) {

	 //        ::numeric::xyzVector<core::Real> _n  = res.xyz("N");
	 //        ::numeric::xyzVector<core::Real> _ca = res.xyz("CA");
	 //        ::numeric::xyzVector<core::Real> _c  = res.xyz("C");

	 //        Eigen::Vector3f n;  n[0]  = _n[0];  n[1]  = _n[1];  n[2] =  _n[2];
	 //        Eigen::Vector3f ca; ca[0] = _ca[0]; ca[1] = _ca[1]; ca[2] = _ca[2];
	 //        Eigen::Vector3f c;  c[0]  = _c[0];  c[1]  = _c[1];  c[2] =  _c[2];

		// 	from_n_ca_c(n,ca,c);
		// }

		PMActor(
			PMActor const & actor0,
			Position const & moveby
		){
			aa_ = actor0.aa_;
			ss_ = actor0.ss_;
			index_ = actor0.index_;			
			set_position(moveby*actor0.position());
		}	

		void 
		set_position(
			Position const & pos
		){ position_ = pos; }

		void 
		moveby(
			Position const & pos
		){ position_ = pos * position_; }

		Position const &
		position() const { return position_; }

		bool operator==(THIS const & o) const {
			return o.position_.translation()==position_.translation() &&
			 	   o.position_.rotation()==position_.rotation() &&	
			       o.aa_      ==aa_       &&
			       o.ss_      ==ss_       &&
			       o.index_   ==index_;
		}

		template<class V>
		void get_n_ca_c_coreReal( V & n, V & ca, V & c ) const {
			V3 tmpn  = position_ * V3( 2.80144, -0.992889, -1.52486 );
		 	V3 tmpca = position_ * V3( 1.95280,  0.220007, -1.52486 );
		 	V3 tmpc  = position_ * V3( 2.87767,  1.4329  , -1.52486 );
		 	n [0] = static_cast<core::Real>(tmpn [0]); n [1] = static_cast<core::Real>(tmpn [1]); n [2] = static_cast<core::Real>(tmpn [2]);
		 	ca[0] = static_cast<core::Real>(tmpca[0]); ca[1] = static_cast<core::Real>(tmpca[1]); ca[2] = static_cast<core::Real>(tmpca[2]);
		 	c [0] = static_cast<core::Real>(tmpc [0]); c [1] = static_cast<core::Real>(tmpc [1]); c [2] = static_cast<core::Real>(tmpc [2]);
		}


		///@brief necessary for testing only
		// bool operator<(THIS const & o) const { return std::make_pair(position_,aa_) < std::make_pair(o.position_,o.aa_); }

		///@brief necessary for testing only
	  	// template<class Archive> void serialize(Archive & ar, const unsigned int ){
	  	// 	ar & position_;
	  	// 	ar & aa_;
	  	// 	ar & ss_;	  		
	  	// 	ar & index_;
	  	// }

	};


template< class P, class MetaData >
void write_pdb( std::ostream & out, PMActor<P> const & a, MetaData const & meta ){
	// typedef Eigen::Matrix<typename P::Scalar,3,1> V3;
	// V3 n,ca,c;
	// a.get_n_ca_c(n,ca,c);
	// int rnum = 1;
	// int anum = 1;
	// chemical::AtomData  ndata( " N  ", "GLY", 'A', rnum, anum+1, "N" );
	// chemical::AtomData cadata( " CA ", "GLY", 'A', rnum, anum+2, "C" );
	// chemical::AtomData  cdata( " C  ", "GLY", 'A', rnum, anum+3, "C" );
	// io::dump_pdb_atom( out,  n,  ndata );
	// io::dump_pdb_atom( out, ca, cadata );
	// io::dump_pdb_atom( out,  c,  cdata );

	// AtomData(
	// 	std::string const & _atomname = AtomData::default_atomname(),        
	// 	std::string const & _resname  = AtomData::default_resname(),       
	// 	char                _chain    = AtomData::default_chain(),     
	// 	int                 _resnum   = AtomData::default_resnum(),      
	// 	int                 _atomnum  = AtomData::default_atomnum(),       
	// 	std::string const & _elem     = AtomData::default_elem(),    
	// 	bool                _ishet    = AtomData::default_ishet(),     
	// 	float               _occ      = AtomData::default_occ(),   
	// 	float               _bfac     = AtomData::default_bfac()
	// );

}



template<class X>
std::ostream & operator<<(std::ostream & out,PMActor<X> const& a){
	return out << "SCPMActor " << a.aa_ << " " << a.index_;
}

}
}

#endif
