#include <QApplication>
#include "mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QFont font = QApplication::font();
    font.setPointSize(10);
    QApplication::setFont(font);

    //QTextCodec *utfCodec = QTextCodec::codecForName("UTF-8"); //creating new utf-8 codec
    //QTextCodec::setCodecForTr(utfCodec);  // setting the utf-8 codec for the tr() tags

    MainWindow w;
    w.showMaximized();
    
    return a.exec();
}
