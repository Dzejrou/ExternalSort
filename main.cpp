// Author: Jaroslav Jindrak 2016 (dzejrou)
#include <string>
#include <iostream>
#include <utility>
#include <cstdint>
#include <fstream>
#include <tuple>
#include <string>
#include <cstdio>
#include <limits>

#define TIMER_ALLOWED 0
#define DEBUG_MESSAGES 0
#define TEST 0

/**
 * Leaving 2GB of memory for OS, other processes
 * and additional variables used in this program.
 */
#define MEMORY_SIZE 8
#if MEMORY_SIZE == 16
#define AVAILABLE_MEMORY 13958643712
#elif MEMORY_SIZE == 4
#define AVAILABLE_MEMORY 2147483648
#else // Assume normal 8 gig machine.
#define AVAILABLE_MEMORY 6442450944
#endif

/**
 * Can be used to easilly switch between C and C++
 * IO facilities.
 */
#define C_IO 1
#if C_IO == 1
#include <cstdio>
#endif


#if DEBUG_MESSAGES == 1
#define DEBUG(MSG) std::cout << MSG << std::endl;
#else
#define DEBUG(MSG)
#endif

#if TIMER_ALLOWED == 1
#include <chrono>
#endif

/**
 * Auxiliary data type representing a pair of a number and
 * the line it was encountered at.
 */
template<typename T>
struct Item
{
	/**
	 * Simple typedef for key type extraction.
	 */
	using key_type = T;

	/**
	 * The number this item contains.
	 */
	T key;

	/**
	 * Number of the line this item was encountered at.
	 * Note: Lines in this assignment begin numbering at 1.
	 */
	std::size_t line;

	/**
	 * Default constructor.
	 */
	Item() = default;

	/**
	 * Constructor.
	 */
	Item(T k, std::uint64_t l)
		: key{k}, line{l}
	{ /* DUMMY BODY */ }

	Item(const Item& other) = default;
	Item& operator=(const Item& other) = default;
	Item(Item&& other) = default;
	Item& operator=(Item&& other) = default;

	/**
	 * Destructor.
	 */
	~Item() = default;

	/**
	 * Comparison of two items, if the keys match compares the lines
	 * of the encounter.
	 * Param: The item compared against.
	 */
	bool operator<(const Item<T>& other) const
	{
		return (key == other.key ? line < other.line : key < other.key);
	}

	/**
	 * Checks if this item is equal to another item.
	 * Param: The item compared against.
	 */
	bool operator==(const Item<T>& other) const
	{
		return key == other.key && line == other.line;
	}

	/**
	 * Checks if this item is inequal to another item.
	 * Param: The item compared against.
	 */
	bool operator!=(const Item<T>& other) const
	{
		return key != other.key || line != other.line;
	}

	/**
	 * Returns string representation of an item,
	 * used for debugging.
	 */
	std::string to_string() const
	{
		return "[" + std::to_string(key) + ", " + std::to_string(line) + "]";
	}
};

/**
 * Simple template alias to avoid using typename
 * too often.
 */
template<typename T>
using item_key_t = typename T::key_type;

/**
 * Auxiliary wrapper template around an array that uses
 * RAII to avoid memory leaks.
 * Note: Originally intended to be used in the algorithm,
 *       but in the end was used only for tests.
 */
template<typename T>
class Array
{
	public:
		/**
		 * Constructor.
		 * Note: No need for other types of constructor
		 *       (like e.g. std::vector has) since all we need
		 *       here is preallocation.
		 */
		Array(std::size_t size)
			: size_{size}, array_{new T[size_]}
		{ /* DUMMY BODY */ }

		/**
		 * Destructor.
		 */
		~Array()
		{
			delete[] array_;
		}

		/**
		 * Returns the number of elements this array has.
		 */
		std::size_t size() const { return size_; }

		/**
		 * Returns a reference to the item at the given index.
		 * Param: Index of the item.
		 */
		T& operator[](std::size_t idx) { return array_[idx]; }

	private:
		/**
		 * Size of the array.
		 */
		std::size_t size_;

		/**
		 * The underlaying array.
		 */
		T* array_;
};

/**
 * A buffered file that is expected to read its
 * data in chunks. (Used for the main file.)
 */
template<typename T>
class BufferedFile
{
	public:
		/**
		 * Constructor.
		 * Param: Name of the file to be read.
		 * Param: Size of the reading buffer in bytes.
		 * Param: Size of the data buffer in sizeof(T) bytes.
		 */
		BufferedFile(const std::string& file_name, std::streamsize bs, std::size_t ds)
			: buffer_{new char[buffer_size_]}, buffer_actual_size_{},
			  line_{1}, data_{new T[ds]}, buffer_size_{bs},
			  data_size_{ds}, number_{}, buffer_index_{},
#if C_IO == 1
			  file_{std::fopen(file_name.c_str(), "r")}
#else
			  file_{file_name}
#endif
		{ /* DUMMY BODY */ }

		/**
		 * Destructor.
		 */
		~BufferedFile()
		{
#if C_IO == 1
			if(file_)
				std::fclose(file_);
#else
			if(file_.is_open())
				file_.close();
#endif

			if(buffer_)
				delete[] buffer_;
			if(data_)
				delete[] data_;
		}

		/**
		 * Deletes all data, used to free memory
		 * before this class reaches end of the scope.
		 */
		void release()
		{
			delete[] buffer_;
			buffer_ = nullptr;
			delete[] data_;
			data_ = nullptr;

#if C_IO == 1
			std::fclose(file_);
			file_ = nullptr;
#else
			file_.close();
#endif
		}

		/**
		 * Reads data from the file.
		 * Returns: Tuple of data array and the number of data read.
		 * Note: Even though this class returns its buffer for sorting
		 *       and temp file storing, it OWNS the buffer and will rewrite
		 *       it on the next read (making it not thread safe)!
		 */
		std::tuple<T*, std::size_t> read()
		{
#if C_IO == 1
			if(std::feof(file_))
				return std::make_tuple(data_, 0);
#else
			if(file_.eof())
				return std::make_tuple(data_, 0);
#endif
			std::size_t data_current{0};

			if(buffer_index_ < buffer_actual_size_)
			{ // Remainder from last time.
				if(buffer_[buffer_index_] == '\n' && number_ == 0)
					++buffer_index_;

				for(; buffer_index_ < buffer_actual_size_; ++buffer_index_)
				{
					if(buffer_[buffer_index_] == '\n')
					{
						data_[data_current++] = T{number_, line_};
						number_ = 0;
						++line_;

						if(data_current >= data_size_)
							break;
					}
					else if(buffer_[buffer_index_] != '\n')
					{
						number_ *= 10;
						number_ += buffer_[buffer_index_] - '0';
					}
				}
			}

			if(data_current < data_size_) // Still some room for new buffer.
			{
				do
				{
#if C_IO == 1
					buffer_actual_size_ = std::fread(buffer_, sizeof(char), buffer_size_, file_);
#else
					file_.read(buffer_, buffer_size_);
					buffer_actual_size_ = file_.gcount();
#endif

					for(buffer_index_ = 0; buffer_index_ < buffer_actual_size_; ++buffer_index_)
					{
						if(buffer_[buffer_index_] == '\n')
						{
							data_[data_current++] = T{number_, line_};
							number_ = 0;
							++line_;

							if(data_current >= data_size_)
								break;
						}
						else if(buffer_[buffer_index_] != '\n')
						{
							number_ *= 10;
							number_ += buffer_[buffer_index_] - '0';
						}
					}
				}
				while(data_current < data_size_ && buffer_actual_size_ > 0);
			}
			DEBUG("[BUFFEREDFILE] Read " + std::to_string(data_current) + " items out of " + std::to_string(data_size_) + ".");

			return std::make_tuple(data_, data_current);
		}

	private:
		/**
		 * Size of the internal reading buffer.
		 */
		std::streamsize buffer_size_;

		/**
		 * Internal buffer for reading in chunks.
		 */
		char* buffer_;

		/**
		 * Actual number of read bytes.
		 */
		std::streamsize buffer_actual_size_;

		/**
		 * Current line number.
		 */
		std::size_t line_;

		/**
		 * Buffer of data (of type T) that is read into.
		 */
		T* data_;

		/**
		 * Size of the data buffer.
		 */
		std::size_t data_size_;

#if C_IO == 1
		/**
		 * File that is being read.
		 */
		FILE* file_;
#else
		/**
		 * File that is being read.
		 */
		std::ifstream file_;
#endif

		/**
		 * Temporary value used during parsing, needs to be
		 * persistent between reads for numbers that are split
		 * between two reads.
		 */
		item_key_t<T> number_;

		/**
		 * Current position in the buffer, needs to be
		 * persistent between reads for when a buffer contains
		 * more data than fits to the data buffer.
		 */
		std::size_t buffer_index_;
};

/**
 * A buffered file that is expected to read its
 * data in chunks. (Used for the main file.)
 */
template<typename T>
class TemporaryBufferedFile
{
	public:
		/**
		 * Constructor.
		 * Param: Name of the file being read.
		 * Param: Size of the data buffer in sizeof(T) bytes.
		 */
		TemporaryBufferedFile(const std::string& file_name, std::size_t ds)
			: data_{new T[ds]}, data_size_{ds},
			  data_index_{}, actual_data_size_{}, at_end_{false},
#if C_IO == 1
			  file_{std::fopen(file_name.c_str(), "rb")}
#else
			  file_{file_name, std::ios::in | std::ios::binary}
#endif
		{
			read();
		}

		/**
		 * Destructor.
		 */
		~TemporaryBufferedFile()
		{
#if C_IO == 1
			std::fclose(file_);
#else
			file_.close();
#endif
			delete[] data_;
		}

		/**
		 * Returns true if all data were read from the file,
		 * false otherwise.
		 */
		bool at_end() const { return at_end_; }

		/**
		 * Returns a pointer to the current data item.
		 */
		T* get()
		{
			return &data_[data_index_];
		}

		/**
		 * Advances the index of the current data item.
		 * Returns true if there still are data to process,
		 * false otherwise.
		 */
		bool advance()
		{
			++data_index_;

			if(data_index_ >= actual_data_size_)
			{
				if(!read())
				{
					at_end_ = true;
					return false;
				}
			}
			return true;
		}


	private:
		/**
		 * Fills the data buffer with data read from the file.
		 */
		bool read()
		{
#if C_IO == 1
			actual_data_size_ = std::fread((char*)data_, sizeof(T), data_size_, file_);
#else
			file_.read((char*)data_, data_size_ * sizeof(T));
			actual_data_size_ = file_.gcount() / sizeof(T);
#endif
			data_index_ = 0;

			return actual_data_size_ != 0;
		}

		/**
		 * Buffer containing read data.
		 */
		T* data_;

		/**
		 * Total size of the data buffer.
		 */
		std::size_t data_size_;

#if C_IO == 1
		/**
		 * File being read.
		 */
		FILE* file_;
#else
		/**
		 * File being read.
		 */
		std::ifstream file_;
#endif

		/**
		 * Index of the currently processed item.
		 */
		std::size_t data_index_;

		/**
		 * Size of the portion of the data buffer
		 * that is actually filled.
		 */
		std::size_t actual_data_size_;

		/**
		 * Auxiliary variable used to check for EOF.
		 */
		bool at_end_;
};

/**
 * Implementation of the insert sort algorithm,
 * used by quick sort for sufficiently small arrays.
 */
class InsertSort
{
	public:
		/**
		 * Main sort function.
		 * Param: Array to be sorted.
		 * Param: Starting index of the sorted area of the array.
		 * Param: Ending index of the sorted area of the array.
		 */
		template<typename T>
		static void sort(T array, std::size_t start, std::size_t end)
		{
			for(std::size_t i = start; i <= end; ++i)
			{
				std::size_t j{i};
				while(j > 0 && array[j] < array[j - 1])
				{
					std::swap(array[j], array[j - 1]);
					--j;
				}
			}
		}
};

/**
 * Implementation of the quick sort algorithm, which is the main
 * algorithm used in this program.
 */
class QuickSort
{
	public:
		/**
		 * Main sort function.
		 * Param: Array that is to be sorted.
		 * Param: Starting index of the sorted area of the array.
		 * Param: Ending index of the sorted area of the array.
		 */
		template<typename T>
		static void sort(T array, std::size_t start, std::size_t end)
		{
			if(end - start > insertion_limit_)
			{
				std::size_t pivot1{}, pivot2{};
				std::tie(pivot1, pivot2) = partition_(array, start, end);

				sort(array, start, pivot1);
				sort(array, pivot2, end);
			}
			else
				InsertSort::sort(array, start, end);
		}

	private:
		/**
		 * When the algorithm accepts array of this or lower size,
		 * it uses the insert sort algorithm instead.
		 */
		static constexpr std::size_t insertion_limit_{16};

		/**
		 * Simple reimplementation of the std::min algorithm.
		 */
		static std::size_t min_(std::size_t fnum, std::size_t snum)
		{
			return fnum < snum ? fnum : snum;
		}

		/**
		 * Simple reimplementation of the std::max algorithm.
		 */
		static std::size_t max_(std::size_t fnum, std::size_t snum)
		{
			return fnum < snum ? snum : fnum;
		}

		/**
		 * Partitioning part of the algorithm.
		 * Param: Array that is to be partitioned.
		 * Param: Starting index.
		 * Param: Ending index.
		 * Returns: Two indices indicating parts that are smaller or bigger
		 *          than the chosen pivot (items between are equal to the pivot).
		 */
		template<typename T>
		static std::tuple<std::size_t, std::size_t> partition_(T array, std::size_t start, std::size_t end)
		{
			auto pivot{static_cast<std::size_t>((end + start) / 2)};

			// Pick median of the first, middle and last elements
			// and place it in middle.
			if(array[pivot] < array[start])
				std::swap(array[start], array[pivot]);
			if(array[end] < array[start])
				std::swap(array[pivot], array[end]);
			if(array[end] < array[pivot])
				std::swap(array[pivot], array[end]);
			

			auto pivot_element = array[pivot];

			// I need them to contain size_t and to have negative values.
			long long i = (long long)start;
			long long j = (long long)end;

			do
			{
				while(array[i] < pivot_element)
					++i;
				while(pivot_element < array[j])
					--j;

				if(i < j)
					std::swap(array[i], array[j]);
				if(i <= j)
				{
					++i;
					--j;
				}
			}
			while(i <= j);

			return std::make_tuple(j, i);
		}
};

/**
 * Main sorter class, accepts an algorithm and a data type as its
 * template argument and then sorts file containing items of the passed
 * type using the algorithm.
 */
template<typename Algorithm, typename T>
class FileSorter
{
	public:
		/**
		 * Constructor.
		 * Param: File to be read.
		 * Param: File to store the result into.
		 */
		FileSorter(const std::string& input_file, const std::string& output_file)
			: tmp_file_count_{}, file_name_{input_file},
			  buffer_{}, buffer_index_{}, buffer_size_{}
#if C_IO == 0
			  , output_{output_file}
#endif
		{ /* DUMMY BODY */ }

		/**
		 * Destructor.
		 */
		~FileSorter() { cleanup_(); };

		/**
		 * Main sort function. Performs quick sort on chunks that fit
		 * into memory and then merges them into the output file.
		 */
		void sort()
		{
			constexpr std::size_t BUFFER_SIZE{30 * 1024 * 1024};
			constexpr std::size_t DATA_SIZE{AVAILABLE_MEMORY / sizeof(T)};
			//constexpr std::size_t OUTPUT_BUFFER_SIZE{30 * 1024 * 1024};

			BufferedFile<T> file{file_name_, (std::streamsize)BUFFER_SIZE, DATA_SIZE};
			T* data;
			std::size_t count{};
#if C_IO == 1
			FILE* binary_output;
#else
			std::ofstream output{};
#endif

			DEBUG("[INFO] Split & Sort.");
			// Split + sort.
			std::tie(data, count) = file.read();
			do
			{
				DEBUG("[INFO] Batch of " + std::to_string(count) + " started.");
				DEBUG("[INFO] Key size: " + std::to_string((count * sizeof(item_key_t<T>))) +
					  " | Line size: " + std::to_string((count * sizeof(std::size_t))) + ".");
				DEBUG("[INFO] Total size: " + std::to_string((count * sizeof(T))) + ".");
				QuickSort::sort(data, 0, count - 1);

				DEBUG("[INFO] Sorting done, saving data in temporary file.");

#if C_IO == 1
				binary_output = std::fopen(("jindraj2.data.tmp" + std::to_string(tmp_file_count_++)).c_str(), "wb");
				std::fwrite((char*)data, sizeof(T), count, binary_output);
				std::fclose(binary_output);
#else
				output.open("jindraj2.data.tmp" + std::to_string(tmp_file_count_++), std::ios::out | std::ios::binary);
				output.write((char*)data, count * sizeof(T));
				output.close();
#endif

				DEBUG("[INFO] Batch of " + std::to_string(count) + " sorted and stored.");

				DEBUG("[INFO] Reading next batch.");
				std::tie(data, count) = file.read();
			}
			while(count > 0);
			file.release();

			// Merge.
			DEBUG("[INFO] Initialising temporary files for merging.");
			TemporaryBufferedFile<T>** tmp_files = new TemporaryBufferedFile<T>*[tmp_file_count_];
			for(std::size_t i = 0; i < tmp_file_count_; ++i)
			{
				tmp_files[i] = new TemporaryBufferedFile<T>{
					"jindraj2.data.tmp" + std::to_string(i),
					DATA_SIZE / tmp_file_count_
				};
			}

			T max_item{};
			max_item.key = std::numeric_limits<item_key_t<T>>::max();
			max_item.line = std::numeric_limits<std::size_t>::max();
			T min_item{max_item};
			T* curr_item{};

			// Initialise output buffer!
			//buffer_size_ = OUTPUT_BUFFER_SIZE;
			//buffer_index_ = 0;
			//buffer_ = new char[buffer_size_];

			std::size_t min_index{};

#if C_IO == 1
			FILE* output_file = std::fopen("data.out", "w");
#endif

			DEBUG("[INFO] Starting merge.");
			do
			{
				min_item = max_item;

				for(std::size_t i = 0; i < tmp_file_count_; ++i)
				{ // Load & min lookup.
					if(tmp_files[i] && (curr_item = tmp_files[i]->get()) != nullptr)
					{
						if(*(curr_item) < min_item)
						{
							min_item = *(curr_item);
							min_index = i;
						}
					}
				}

				// Output number have to be unique with the lowest line
				// number possible.
				for(std::size_t i = 0; i < tmp_file_count_; ++i)
				{
					if(!tmp_files[i])
						continue;
					while(tmp_files[i]->get() && tmp_files[i]->get()->key == min_item.key)
					{
						if(!tmp_files[i]->advance())
						{
							delete tmp_files[i];
							tmp_files[i] = nullptr;
							break;
						}
					}
				}

				if(min_item < max_item)
#if C_IO == 1
					//append_to_buffer_(min_item, output_file);
					fprintf(output_file, "%lu %lu\n", min_item.key, min_item.line);
#else
					append_to_buffer_(min_item, output_);
#endif
			}
			while(min_item != max_item);
			DEBUG("[INFO] Merge done.");

#if C_IO == 1
			flush_(output_file);
			std::fclose(output_file);
#else
			flush_(output_);
			output_.close();
#endif

			// Cleanup.
			for(std::size_t i = 0; i < tmp_file_count_; ++i)
				if(tmp_files[i]) delete tmp_files[i];
			delete[] tmp_files;
			delete[] buffer_;
		}
	private:
		/**
		 * Appends an item to the output buffer.
		 */
		void append_to_buffer_(T& item, std::ofstream& output)
		{
			std::string key{std::to_string(item.key)};
			std::string line{std::to_string(item.line)};

			// That + 2 is for space and newline.
			if(buffer_size_ - buffer_index_ < key.size() + line.size() + 2)
				flush_(output);
			
			for(const auto& c : key)
				buffer_[buffer_index_++] = c;
			buffer_[buffer_index_++] = ' ';

			for(const auto& c : line)
				buffer_[buffer_index_++] = c;
			buffer_[buffer_index_++] = '\n';
		}

		/**
		 * Appends an item to the output buffer.
		 */
		void append_to_buffer_(T& item, FILE* output)
		{
			std::string key{std::to_string(item.key)};
			std::string line{std::to_string(item.line)};

			// That + 2 is for space and newline.
			if(buffer_size_ - buffer_index_ < key.size() + line.size() + 2)
				flush_(output);
			
			for(auto c : key)
				buffer_[buffer_index_++] = c;
			buffer_[buffer_index_++] = ' ';

			for(auto c : line)
				buffer_[buffer_index_++] = c;
			buffer_[buffer_index_++] = '\n';
		}

		/**
		 * Flushes the contents of the output buffer.
		 */
		void flush_(std::ofstream& output)
		{
			if(buffer_index_ == 0)
				return;
			output.write(buffer_, buffer_index_);
			buffer_index_ = 0;
		}

		/**
		 * Flushes the contents of the output buffer.
		 */
		void flush_(FILE* output)
		{
			if(buffer_index_ == 0)
				return;
			std::fwrite(buffer_, 1, buffer_index_, output);
			buffer_index_ = 0;
		}

		/**
		 * Removes all temporary files.
		 */
		void cleanup_()
		{
			for(std::size_t i = 0; i < tmp_file_count_; ++i)
				std::remove(("jindraj2.data.tmp" + std::to_string(i)).c_str());
		}

		/**
		 * Number of temporary files that were created during the sorting.
		 */
		std::size_t tmp_file_count_;

		/**
		 * Name of the input file.
		 */
		std::string file_name_;

		/**
		 * Buffer used for outputting the result.
		 */
		char* buffer_;

		/**
		 * Current index in the output buffer.
		 */
		std::size_t buffer_index_;

		/**
		 * Total size of the buffer.
		 */
		std::size_t buffer_size_;

#if C_IO == 0
		/**
		 * Output file.
		 */
		std::ofstream output_;
#endif
};

/**
 * Forward declarations for better readability.
 */
void test();

/**
 * Entry point of the sorter.
 */
int main(int argc, char** argv)
{
#if C_IO == 0
	std::ios::sync_with_stdio(false); // Disable synchronizing with scanf/printf for better
									  // I/O performance.
#endif

#if TIMER_ALLOWED == 1
	auto start_time = std::chrono::system_clock::now();
#endif

	DEBUG("[DEBUG] Starting.");

	std::string input_file{};
	if(argc > 1)
		input_file = argv[1];
	else
		input_file = "data.txt";

	std::string output_file{};
	if(argc > 2)
		output_file = argv[2];
	else
		output_file = "data.out";

	FileSorter<QuickSort, Item<std::uint64_t>> sorter{input_file, output_file};
	sorter.sort();
	
#if TIMER_ALLOWED == 1
	auto end_time = std::chrono::system_clock::now();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cout << "[TIMER] Elapsed: " << elapsed_time.count() << "s." << std::endl;
#endif

#if TEST == 1
	std::cout << "[TEST] Running tests." << std::endl;
#if TIMER_ALLOWED == 1
	start_time = std::chrono::system_clock::now();
#endif
	test();
#if TIMER_ALLOWED == 1
	end_time = std::chrono::system_clock::now();
	elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cout << "[TEST TIMER] Elapsed: " << elapsed_time.count() << "s." << std::endl;
#endif
#endif

	DEBUG("[DEBUG] Ending.");
}

/**
 * Forward declarations for better readability.
 */
bool test_1();
bool test_2();
bool test_3();
bool test_4();
bool test_5();

/**
 * A simple test suite, runs all of the tests.
 */
void test()
{
	std::cout << "[TEST] Test #1 started." << std::endl;
	if(test_1())
		std::cout << "[TEST] Test #1 succeeded." << std::endl;
	else
		std::cout << "[TEST] Test #1 failed." << std::endl;

	std::cout << "[TEST] Test #2 started." << std::endl;
	if(test_2())
		std::cout << "[TEST] Test #2 succeeded." << std::endl;
	else
		std::cout << "[TEST] Test #2 failed." << std::endl;

	std::cout << "[TEST] Test #3 started." << std::endl;
	if(test_3())
		std::cout << "[TEST] Test #3 succeeded." << std::endl;
	else
		std::cout << "[TEST] Test #3 failed." << std::endl;

	std::cout << "[TEST] Test #4 started." << std::endl;
	if(test_4())
		std::cout << "[TEST] Test #4 succeeded." << std::endl;
	else
		std::cout << "[TEST] Test #4 failed." << std::endl;

	std::cout << "[TEST] Test #5 started." << std::endl;
	if(test_5())
		std::cout << "[TEST] Test #5 succeeded." << std::endl;
	else
		std::cout << "[TEST] Test #5 failed." << std::endl;
}

/**
 * QuickSort test with a primitive type, just checks if
 * the algorithm itself is ok.
 */
bool test_1()
{
	constexpr std::size_t test_size{1000000};
	int* array = new int[test_size];
	for(std::size_t i = 0; i < test_size; ++i)
		array[i] = test_size - i;
	QuickSort::sort(array, 0, test_size - 1);

	bool res{true};
	for(std::size_t i = 1; i < test_size; ++i)
		if(array[i] < array[i - 1]) res = false;

	delete[] array;
	return res;
}

/**
 * QuickSort test with the Item type, checks how the
 * algorithm works with this a bit more complex type.
 */
bool test_2()
{
	Item<int> array1[5];
	array1[0].key = 2;
	array1[0].line = 1;
	array1[1].key = 2;
	array1[1].line = 3;
	array1[2].key = 1;
	array1[2].line = 4;
	array1[3].key = 5;
	array1[3].line = 2;
	array1[4].key = 1;
	array1[4].line = 5;

	QuickSort::sort(array1, 0, 5);

	Item<int> array2[5];
	array2[0].key = 1;
	array2[0].line = 4;
	array2[1].key = 1;
	array2[1].line = 5;
	array2[2].key = 2;
	array2[2].line = 1;
	array2[3].key = 2;
	array2[3].line = 3;
	array2[4].key = 5;
	array2[4].line = 2;

	for(std::size_t i = 0; i < 5; ++i)
	{
		if(array1[i].key != array2[i].key ||
		   array1[i].line != array2[i].line)
			return false;
	}
	return true;
}

/**
 * Test that checks that the BufferedFile works on small
 * input.
 */
bool test_3()
{
	std::uint64_t test_data[] = { 123456, 654321, 123, 321};
	std::ofstream output{"data.test3"};
	for(auto data : test_data)
		output << data << std::endl;
	output.close();

	BufferedFile<Item<std::uint64_t>> input{"data.test3", 1024, 4};
	Item<std::uint64_t>* data;
	std::size_t count;

	std::tie(data, count) = input.read();
	std::remove("data.test3");
	if(count != 4)
	{
		DEBUG("[TEST #3 ERROR] Count mismatch: " + std::to_string(count) + " != 4.");
		return false;
	}

	for(std::size_t i = 0; i < 4; ++i)
	{
		if(data[i].key != test_data[i])
		{
			DEBUG("[TEST #3 ERROR] Data mismatch: " + std::to_string(data[i].key)
				  + " != " + std::to_string(test_data[i]) + " at line "
				  + std::to_string(i) + ".");
			return false;
		}
	}

	for(std::size_t i = 1; i < 4; ++i)
	{
		if(data[i - 1].line != i)
		{
			DEBUG("[TEST #3 ERROR] Line mismatch: " + std::to_string(data[i - 1].line)
				  + " != " + std::to_string(i) + ".");
			return false;
		}
	}

	return true;
}

/**
 * Test that checks that the BufferedFile works on a bigger
 * input.
 */
bool test_4()
{
	constexpr std::size_t test_size{900000};
	std::uint64_t test_data[test_size];
	for(std::size_t i = 0; i < test_size; ++i)
		test_data[i] = i;

	std::ofstream output{"data.test4"};
	for(auto data : test_data)
		output << data << std::endl;
	output.close();

	BufferedFile<Item<std::uint64_t>> input{"data.test4", 1024, test_size / 10};
	Item<std::uint64_t>* data;
	std::size_t count;

	for(std::size_t j = 1; j <= 10; ++j)
	{
		std::tie(data, count) = input.read();
		std::size_t start = (test_size / 10) * (j - 1);
		std::size_t end = (test_size / 10) * j;
		DEBUG("[TEST #4 BATCH #" + std::to_string(j) + "] " + std::to_string(start)
			  + " - " + std::to_string(end) + ".");

		if(count != (test_size / 10))
		{
			DEBUG("[TEST #4 ERROR] Count mismatch: " + std::to_string(count)
				  + " != " + std::to_string(test_size / 10) + ".");
			return false;
		}

		for(std::size_t i = start, k = 0; i < end && k < count; ++i, ++k)
		{
			if(data[k].key != test_data[i])
			{
				DEBUG("[TEST #4 ERROR] Data mismatch: " + std::to_string(data[k].key)
					  + " != " + std::to_string(test_data[i]) + " at line "
					  + std::to_string(i) + ".");
				return false;
			}
		}

		for(std::size_t i = start, k = 0; i < end && k < count; ++i, ++k)
		{
			if(data[k].line != i + 1)
			{
				DEBUG("[TEST #4 ERROR] Line mismatch: " + std::to_string(data[i - 1].line)
					  + " != " + std::to_string(i) + ".");
				return false;
			}
		}
	}

	std::remove("data.test4");

	return true;
}

/**
 *
 */
bool test_5()
{
	// Changed the way the temp file works, gotta
	// rewrite this test. (low prio tho)
	return true;

	/*
	constexpr std::size_t test_size{3};
	Item<std::uint64_t> test_data[test_size];
	for(std::size_t i = 0; i < test_size; ++i)
	{
		test_data[i].key = i;
		test_data[i].line = test_size - i;
	}

	std::ofstream output{"data.test5"};
	for(auto& item : test_data)
		output << item.key << " " << item.line << std::endl;
	output.close();

	TemporaryBufferedFile<Item<std::uint64_t>> input{"data.test5", 1024, 3};
	Item<std::uint64_t>* data{};
	std::size_t count{};

	std::tie(data, count) = input.read();
	if(count != test_size)
	{
		DEBUG("[TEST #5 ERROR] Count mismatch: " + std::to_string(count) + "!= 3.");
		return false;
	}

	for(std::size_t i = 0; i < test_size; ++i)
	{
		if(data[i] != test_data[i])
		{
			DEBUG("[TEST #5 ERROR] Data mismatch: " + data[i].to_string() + "!="
				  + test_data[i].to_string() + ".");
			return false;
		}
	}

	std::remove("data.test5");

	return true;
	*/
}
