// RaMA integration support. Serialization declarations retained for RMQ compatibility.
#pragma once
#include "gsacak.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
namespace ramag::rama_source {
inline bool canonical(unsigned char c) { return c=='A'||c=='C'||c=='G'||c=='T'; }
class Serializable {
public:
	virtual void serialize(std::ostream& out) const = 0;
	virtual void deserialize(std::istream& in) = 0;

	bool saveToFile(const std::string& filename) const;
	bool loadFromFile(const std::string& filename);
protected:
	template<typename T>
	void saveNumber(std::ostream& out, const T& value) const {
		out.write(reinterpret_cast<const char*>(&value), sizeof(T));
	}

	template<typename T>
	void loadNumber(std::istream& in, T& value) {
		in.read(reinterpret_cast<char*>(&value), sizeof(T));
	}

	template<typename T>
	void saveArray(std::ostream& out, const T* array, size_t size) const {
		out.write(reinterpret_cast<const char*>(array), sizeof(T) * size);
	}

	template<typename T>
	void loadArray(std::istream& in, T* array, size_t size) {
		in.read(reinterpret_cast<char*>(array), sizeof(T) * size);
	}

	template<typename T>
	void saveVector(std::ostream& out, const std::vector<T>& vec) const {
		size_t size = vec.size();
		saveNumber(out, size);
		for (const T& item : vec) {
			saveNumber(out, item);
		}
	}

	template<typename T>
	void loadVector(std::istream& in, std::vector<T>& vec) {
		size_t size;
		loadNumber(in, size);
		vec.resize(size);
		for (T& item : vec) {
			loadNumber(in, item);
		}
	}

	template<typename T>
	void saveVector2D(std::ostream& out, const std::vector<std::vector<T>>& vec2d) const {
		size_t size = vec2d.size();
		saveNumber(out, size);
		for (const auto& vec : vec2d) {
			saveVector(out, vec);
		}
	}

	template<typename T>
	void loadVector2D(std::istream& in, std::vector<std::vector<T>>& vec2d) {
		size_t size;
		loadNumber(in, size);
		vec2d.resize(size);
		for (auto& vec : vec2d) {
			loadVector(in, vec);
		}
	}
};


struct SequenceInfo {
    std::string_view sequence;
    std::string_view header;
    uint_t seq_len;
};
template<class T> T getMinValue(const T& a,const T& b) { return a<b?a:b; }
template<class T> T getMaxValue(const T& a,const T& b) { return a>b?a:b; }
inline std::string joinPaths(const std::string& a,const std::string& b) {
    return (std::filesystem::path(a)/b).string();
}
}

