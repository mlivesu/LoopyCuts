#ifndef PATH_ANISOTROPIC_GEODESIC
#define PATH_ANISOTROPIC_GEODESIC

    /**
     * @brief The Path struct stores information of a path and its total length.
     */
    template <class ScalarType>
    struct Path {
        ScalarType            length;   ///< The length of the whole path.
        std::vector<size_t>   nodes;    ///< The sequence of nodes of the path from the first to the last.
        bool isTJunction;

        /**
       * @brief Path Default constructor.
       */
        Path() : length(std::numeric_limits<ScalarType>::max()) { isTJunction=false;}

        /**
       * @brief Path Constructor initializer.
       * @param _length
       */
        Path(const ScalarType _length) : length(_length) {isTJunction=false; }

        /**
       * @brief Path Copy constructor.
       * @param p
       */
        Path(const Path &p) : length(p.length), nodes(p.nodes),isTJunction(p.isTJunction) { }

        /**
       * @brief operator = Copy assignment.
       * @param p
       * @return
       */
        Path & operator =(const Path &p) {
            if (this != &p) {
                length = p.length;
                nodes = p.nodes;
                isTJunction=p.isTJunction;
            }
            return *this;
        }

#if __cplusplus >= 201103L
        /**
       * @brief Path Move constructor.
       * @param p
       */
        Path(Path &&p) : length(p.length),
            nodes(std::move(p.nodes)),
            isTJunction(p.isTJunction)
        {
            p.length = std::numeric_limits<ScalarType>::max();
            p.nodes.clear();
            p.isTJunction=false;
        }

        /**
       * @brief operator = Move assignment.
       * @param p
       * @return
       */
        Path & operator =(Path &&p) {
            if (this != &p) {
                length = p.length;
                nodes = std::move(p.nodes);
                p.length = std::numeric_limits<ScalarType>::max();
                p.nodes.clear();
                p.isTJunction=false;
            }
            return *this;
        }
#endif

        /**
       * @brief clear
       */
        void clear() {
            length = std::numeric_limits<ScalarType>::max();
            nodes.clear();
        }

        /**
       * @brief operator <
       * @param p
       * @return
       */
        bool operator <(const Path &p) const {
            if (length != p.length)
                return length < p.length;
            if (nodes.size() != p.nodes.size())
                return nodes.size() < p.nodes.size();
            return nodes < p.nodes;
        }

        bool operator >(const Path &p) const {
            return !(*this < p);
        }
        bool operator <=(const Path &p) const {
            return !(p < *this);
        }
        bool operator >=(const Path &p) const {
            return !(*this < p);
        }
        bool operator ==(const Path &p) const {
            return !(*this < p) && !(p < *this);
        }
        bool operator !=(const Path &p) const {
            return *this < p || p < *this;
        }

    };

#endif //PATH_ANISOTROPIC_GEODESIC
