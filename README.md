# 🎧 Dijital Ses İşleme Stüdyosu v2.0

Gelişmiş ses düzenleme, kullanıcı yönetimi, FFT ile spektrum analizi, MIDI klavye simülasyonu ve şifreli ses kaydı özelliklerine sahip C++ ile geliştirilmiş bir masaüstü uygulaması.

---

## 🇹🇷 Türkçe Açıklama

**Dijital Ses İşleme Stüdyosu v2.0**, WAV ve MP3 dosyaları üzerinde çeşitli ses işlemleri ve analizleri yapabilen bir yazılımdır.

### ✨ Özellikler

- 🔐 Kullanıcı yönetimi ve şifreleme (XOR)
- 🎵 WAV ve MP3 dosyası yükleme/kaydetme
- 🎚️ Echo, Reverb, Distortion, Normalize, Fade efektleri
- 📊 FFT tabanlı Spektrum Analizi
- 🔁 Çoklu Kanal Karıştırma
- 🎹 MIDI Klavye Simülasyonu (Beep ile nota çalımı)
- 🔐 Şifreli ses verisi kaydetme ve yükleme

---

## 🇬🇧 English Description

**Digital Audio Processing Studio v2.0** is a C++ based desktop application for WAV/MP3 audio editing, encryption, spectrum analysis, and MIDI simulation.

### ✨ Features

- 🔐 Encrypted user login system (XOR encryption)
- 🎵 Load and save WAV/MP3 files
- 🎚️ Apply Echo, Reverb, Distortion, Normalize, and Fade
- 📊 Real-time FFT Spectrum Analysis
- 🔁 Multi-channel mixing with pan and volume controls
- 🎹 Simulated MIDI keyboard using PC keys and Windows Beep
- 🔐 Secure audio file save/load with encryption

---

## 📁 Gereksinimler / Requirements

- Windows OS
- C++11 destekli bir derleyici (MinGW, MSVC)
- `dr_wav.h` ve `dr_mp3.h` başlık dosyaları (dahil)
- `conio.h`, `windows.h` (Windows spesifik)

---

## 🚀 Derleme Talimatları / Build Instructions

```bash
g++ main.cpp -o AudioStudio -std=c++11
```

> Not: `dr_wav.h` ve `dr_mp3.h` dosyaları ile aynı dizinde olduğundan emin olun.

---

## 🧠 Eğitim Amaçlı Kullanım

Bu proje ses işleme algoritmaları, kullanıcı yönetimi, şifreleme ve FFT gibi konuları öğrenmek isteyen öğrenciler için uygundur.