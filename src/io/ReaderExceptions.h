//
// Created by jan on 15/10/19.
//

#ifndef LFNS_READEREXCEPTION_H
#define LFNS_READEREXCEPTION_H

namespace io {
    class ReaderException : public std::exception {
    public:
        ReaderException(std::string error_message) : _error_message(error_message) {}

        const char *what() const throw() { return _error_message.c_str(); }

    private:
        std::string _error_message;
    };

    class EntryNotFoundException : public ReaderException {
    public:
        EntryNotFoundException(std::string error_message) : ReaderException(error_message) {}
    };
}

#endif //LFNS_READEREXCEPTION_H
