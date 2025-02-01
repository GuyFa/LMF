/**
 * Simple formatted-logger class with built-in support for indentation and
 * outputting to a user-specified stream.
 *
 * To see an example of how it works, check out test/test_logger.cpp.
 */

#pragma once

#include <string>
#include <fstream>
#include <memory>
#include <sstream>

class Logger : public std::ostream
{
public:
    Logger();
    Logger(int log_level);
    ~Logger();

	int& indent() {
		return m_indent;
	}
    int level() const;
    int set_level(int log_level);
	void restore_cout();
	bool quiet() const;
	bool quiet(bool q);
    std::string &set_normal_prefix(const char *log_prefix);
    const std::string &normal_prefix() const;
    std::string &set_warn_prefix(const char *log_prefix);
    const std::string &warn_prefix() const;
    std::string &set_error_prefix(const char *log_prefix);
    const std::string &error_prefix() const;
    bool auto_flush() const;
    bool set_auto_flush(bool auto_flush);
    Logger& flush();

	// http://forums.cgsociety.org/showthread.php?f=89&t=1322430
	void set_send_to_viewer(int bSend) {
		m_bSendToViewer = bSend;
	}
	void send_to_viewer(std::string msg, bool bAddSeparator);
	void send_color_to_viewer(std::string color);

    // returns 0 on success
    int use_file(const char *fname, bool append = false);
    void use_stream(std::ostream &out);

    // formatted output versions
    Logger& start_block(int level, const char *format, ...);
    Logger& end_block(int level, const char *format, ...);
    Logger& print(const char *format, ...);
	Logger& println(const char *format, ...);
    Logger& println(int level, const char *format, ...);
    Logger& start_block(const char *format, ...);
    Logger& end_block(const char *format, ...);
	Logger& emphasize(const char *format, ...);
	Logger& emphasize2(const char *format, ...);
    Logger& warn(const char *format, ...);
    Logger& error(const char *format, ...);

    // ostream versions
    std::ostream& print(int level);
    std::ostream& print();
    std::ostream& start_block(int level);
    std::ostream& start_block();
    std::ostream& end_block(int level);
    std::ostream& end_block();
	std::ostream& emphasize();
	std::ostream& emphasize2();
    std::ostream& warn();
    std::ostream& error();

	template<typename T>
    Logger& operator<<(const T &v) {
		std::stringstream ss;
		ss << v;
		(*m_out) << ss.str();
		send_to_viewer(ss.str(), 1);
        return *this;
    }
    // for std::endl and other manipulators
    typedef std::ostream& (*STRFUNC)(std::ostream&);
    Logger& operator<<(STRFUNC f) {
		std::stringstream ss;
		f(ss);
		(*m_out) << ss.str();
		send_to_viewer(ss.str(), 1);
        return *this;
    }

    // This class is useful for auto-indenting within code blocks and
    // methods.
    class Block {
    protected:
        Block(Logger &log, int log_level, const char *message);
        //Block(Block &block);
        Logger *m_pLog;
        bool m_enabled;
        char m_message[1025];

        friend class Logger;
    public:
        ~Block();
    };

    // create a block object
    Block block(int log_level, const char *message);
    Block block(int log_level, const std::string &message);

protected:
    enum { MSG_NORMAL = 0, MSG_WARN = 1, MSG_ERROR = 2 };

    int m_indent;  // indentation level
    int m_logging_level;  // logging level
	bool m_quiet; // be quiet
	std::streambuf *m_cout_sbuf, *m_cerr_sbuf;
    std::string m_normal_prefix; // prefix (displayed before indentation)
    std::string m_warn_prefix; // prefix (displayed before indentation)
    std::string m_error_prefix; // prefix (displayed before indentation)
    std::ostream *m_out;  // logging destination
	std::stringstream m_stbuf;
    std::auto_ptr<std::ostream> m_null; // null IO stream (useful for discarding)
    bool m_clean_up_file; // whether the file needs to be closed and delete-ed
    bool m_auto_flush;    // whether to flush the output every time
	int m_bSendToViewer; // 0 - no; 1 - yes, but change to zero if a problem occurs; 2 - yes

    void initialize();
    void deinitialize();
    void clean_up_file();
    void output_header();
    void vprintf_helper(const char *message_prefix, int type, const char *format, va_list ap, bool bNewLine = true);
    void output_prefix(const char *message_prefix, int type);

private:
    // disallow copying
    Logger(const Logger &log);
    Logger &operator=(const Logger &log);

	// for ostream capabilities
	struct CLogBuf : public std::stringbuf {
        CLogBuf(Logger &log) : m_log(log) {}
        ~CLogBuf() {  pubsync(); }
        int sync() {
            (*m_log.m_out) << str();
			m_log.send_to_viewer(str(), 1);
            str("");
            return (*m_log.m_out) ? 0 : -1;
        }

	private:
		Logger &m_log;
	};
};



