#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <conio.h>
#include <windows.h>
#include <complex>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"
#define DR_MP3_IMPLEMENTATION
#include "dr_mp3.h" 

using namespace std;

struct User {
    string username;
    string password;
    bool isAdmin;
};

vector<User> users;
User currentUser;
bool loggedIn = false;

void saveUsersToFile(); 


struct AudioData {
    vector<float> samples;
    int sampleRate = 44100;
    int channels = 1;
    int totalSamples = 0;
    string filename = "";
};
AudioData audio;
bool audioLoaded = false;

// Gelismis sifreleme (basit XOR sifreleme)
const string ENCRYPTION_KEY = "SesStudyosuGizliAnahtar123";

string encrypt(const string& input) {
    string out = input;
    for (size_t i = 0; i < out.length(); i++) {
        out[i] = out[i] ^ ENCRYPTION_KEY[i % ENCRYPTION_KEY.length()];
    }
    return out;
}
string decrypt(const string& input) {
    return encrypt(input); // XOR simetrik
}

// FFT Uygulamasi
void fft(vector<complex<double>>& x) {
    const size_t N = x.size();
    if (N <= 1) return;

    // Cift ve tek indisli elemanları ayarı
    vector<complex<double>> even(N/2), odd(N/2);
    for (size_t i = 0; i < N/2; i++) {
        even[i] = x[2*i];
        odd[i] = x[2*i+1];
    }

    // Tekrarla FFT uygula
    fft(even);
    fft(odd);

    // Birlestir
    for (size_t k = 0; k < N/2; k++) {
        complex<double> t = exp(complex<double>(0, -2 * M_PI * k / N)) * odd[k];
        x[k] = even[k] + t;
        x[k + N/2] = even[k] - t;
    }
}

bool login() {
    string u, p;
    int attempts = 0;

    while (attempts < 3) {
        cout << "Kullanici Adi: ";
        cin >> u;
        cout << "Sifre: ";

        p = "";
        char ch;
        while ((ch = _getch()) != 13) {
            if (ch == 8 && !p.empty()) {
                p.pop_back();
                cout << "\b \b";
            } else {
                p += ch;
                cout << "*";
            }
        }
        cout << endl;

        for (const auto& user : users) {
            if (user.username == u && user.password == encrypt(p)) {
                currentUser = user;
                loggedIn = true;
                cout << "Hos geldiniz, " << currentUser.username << " (" << (currentUser.isAdmin ? "Admin" : "User") << ")\n";
                return true;
            }
        }

        attempts++;
        cout << "Hatali giris. Kalan hak: " << 3 - attempts << "\n";
    }

    cout << "3 defa hatali giris. Programdan cikiliyor.\n";
    return false;
}

// WAV dosyası yükleme fonksiyonu
bool loadWavFile() {
    string filename;
    cout << "WAV dosya adini girin (orn: ses.wav): ";
    cin >> filename;
    
    // dr_wav ile wav dosyasını yükleme
    drwav wav;
    if (!drwav_init_file(&wav, filename.c_str(), NULL)) {
        cout << "[-] Hata: WAV dosyasi acilamadi.\n";
        return false;
    }
    
    // Dosya bilgilerini kaydet
    audio.sampleRate = wav.sampleRate;
    audio.channels = wav.channels;
    audio.totalSamples = wav.totalPCMFrameCount;
    audio.filename = filename;
    
    // Örnekleri oku
    audio.samples.resize(wav.totalPCMFrameCount * wav.channels);
    drwav_read_pcm_frames_f32(&wav, wav.totalPCMFrameCount, audio.samples.data());
    
    // Kaynakları temizle
    drwav_uninit(&wav);
    
    // Ses yüklendi
    audioLoaded = true;
    
    // Bilgi gÃ¶ster
    float duration = (float)audio.totalSamples / audio.sampleRate;
    cout << "[+] WAV dosyasi yuklendi: " << filename << "\n";
    cout << "    Ornekleme hizi: " << audio.sampleRate << " Hz\n";
    cout << "    Kanal sayisi: " << audio.channels << "\n";
    cout << "    Toplam ornek: " << audio.totalSamples << "\n";
    cout << "    Sure: " << duration << " saniye\n";
    
    return true;
}

// MP3 dosyasını yükleme fonksiyonu 
bool loadMp3File() {
    string filename;
    cout << "MP3 dosya adini girin (orn: muzik.mp3): ";
    cin >> filename;
    
    // dr_mp3 ile MP3 dosyasını yükleme
    drmp3 mp3;
    if (!drmp3_init_file(&mp3, filename.c_str(), NULL)) {
        cout << "[-] Hata: MP3 dosyasi acilamadi.\n";
        return false;
    }
    
    // MP3 bilgilerini al
    drmp3_uint64 totalFrames = drmp3_get_pcm_frame_count(&mp3);
    
    // Dosya bilgilerini kaydet
    audio.sampleRate = mp3.sampleRate;
    audio.channels = mp3.channels;
    audio.totalSamples = totalFrames;
    audio.filename = filename;
    
    // Örnekleri oku
    audio.samples.resize(totalFrames * mp3.channels);
    drmp3_read_pcm_frames_f32(&mp3, totalFrames, audio.samples.data());
    
    // Kaynakları temizle
    drmp3_uninit(&mp3);
    
    // Ses yüklendi
    audioLoaded = true;
    
    // Bilgi göster
    float duration = (float)audio.totalSamples / audio.sampleRate;
    cout << "[+] MP3 dosyasi yuklendi: " << filename << "\n";
    cout << "    Ornekleme hizi: " << audio.sampleRate << " Hz\n";
    cout << "    Kanal sayisi: " << audio.channels << "\n";
    cout << "    Toplam ornek: " << audio.totalSamples << "\n";
    cout << "    Sure: " << duration << " saniye\n";
    
    return true;
}

// WAV dosyası kaydetme
bool saveWavFile() {
    if (!audioLoaded) {
        cout << "[-] Kaydedilecek ses verisi bulunamadi.\n";
        return false;
    }
    
    string filename;
    cout << "Kaydedilecek WAV dosya adini girin (orn: kayit.wav): ";
    cin >> filename;
    
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = audio.channels;
    format.sampleRate = audio.sampleRate;
    format.bitsPerSample = 32;
    
    drwav wav;
    if (!drwav_init_file_write(&wav, filename.c_str(), &format, NULL)) {
        cout << "[-] Hata: WAV dosyasi olusturulamadi.\n";
        return false;
    }
    
    drwav_uint64 framesWritten = drwav_write_pcm_frames(&wav, audio.totalSamples, audio.samples.data());
    drwav_uninit(&wav);
    
    if (framesWritten != audio.totalSamples) {
        cout << "[-] Uyari: Bazi ornekler yazilamadi.\n";
    }
    
    cout << "[+] WAV dosyasi basariyla kaydedildi: " << filename << "\n";
    return true;
}

// Sahte ses uretme
void generateDummyAudio() {
    int sure;
    int frekans;
    cout << "Ses suresi (saniye): ";
    cin >> sure;
    cout << "Frekans (Hz): ";
    cin >> frekans;
    
    audio.samples.clear();
    audio.sampleRate = 44100;
    
    for (int i = 0; i < audio.sampleRate * sure; i++) {
        float val = sin(2 * M_PI * frekans * i / audio.sampleRate);
        audio.samples.push_back(val);
    }
    audio.channels = 1;
    audio.totalSamples = audio.sampleRate * sure;
    audio.filename = "generated.wav";
    audioLoaded = true;
    cout << "[+] Ses simulasyonu olusturuldu (" << sure << " saniye, " << frekans << "Hz).\n";
}

// Gelismis Sifreli WAV veri kaydetme
void saveAudio() {
    if (!audioLoaded) {
        cout << "[-] Kaydedilecek ses verisi bulunamadi.\n";
        return;
    }
    
    string filename;
    cout << "Sifreli ses dosyasi adini girin: ";
    cin >> filename;
    
    ofstream out(filename, ios::binary);
    
    // Metaveri yazma (sample rate, kanal sayisi, toplam ornek)
    out.write(reinterpret_cast<char*>(&audio.sampleRate), sizeof(int));
    out.write(reinterpret_cast<char*>(&audio.channels), sizeof(int));
    out.write(reinterpret_cast<char*>(&audio.totalSamples), sizeof(int));
    
    // Ornekleri sifrele ve yaz
    for (float sample : audio.samples) {
        // XOR shifreleme icin anahtari donguluyoruz
        for (size_t i = 0; i < sizeof(float); i++) {
            reinterpret_cast<char*>(&sample)[i] ^= ENCRYPTION_KEY[i % ENCRYPTION_KEY.length()];
        }
        out.write(reinterpret_cast<char*>(&sample), sizeof(float));
    }
    
    out.close();
    cout << "[+] Ses verisi dosyaya sifreli olarak kaydedildi: " << filename << "\n";
}

// Gelismis Sifreli WAV veri yukleme
void loadAudio() {
    string filename;
    cout << "Sifreli ses dosyasi adini girin: ";
    cin >> filename;
    
    ifstream in(filename, ios::binary);
    if (!in) {
        cout << "[-] Dosya acilamadi veya bulunamadi.\n";
        return;
    }
    
    // Metaveri okuma
    in.read(reinterpret_cast<char*>(&audio.sampleRate), sizeof(int));
    in.read(reinterpret_cast<char*>(&audio.channels), sizeof(int));
    in.read(reinterpret_cast<char*>(&audio.totalSamples), sizeof(int));
    
    // Ornekleri oku ve sifresini coz
    audio.samples.clear();
    float sample;
    while (in.read(reinterpret_cast<char*>(&sample), sizeof(float))) {
        // XOR sifre cozme
        for (size_t i = 0; i < sizeof(float); i++) {
            reinterpret_cast<char*>(&sample)[i] ^= ENCRYPTION_KEY[i % ENCRYPTION_KEY.length()];
        }
        audio.samples.push_back(sample);
    }
    
    in.close();
    audioLoaded = true;
    audio.filename = filename;
    cout << "[+] Sifreli ses verisi basariyla yuklendi (" << audio.samples.size() << " ornek).\n";
    cout << "    Ornekleme hizi: " << audio.sampleRate << " Hz\n";
    cout << "    Kanal sayisi: " << audio.channels << "\n";
    cout << "    Toplam ornek: " << audio.totalSamples << "\n";
    cout << "    Sure: " << (float)audio.totalSamples / audio.sampleRate << " saniye\n";
}

// Ses bilgilerini görüntüle
void displayAudioInfo() {
    if (!audioLoaded) {
        cout << "[-] Yuklenmis ses dosyasi bulunamadi.\n";
        return;
    }
    
    float duration = (float)audio.totalSamples / audio.sampleRate;
    cout << "\n=== Ses Dosyasi Bilgileri ===\n";
    cout << "Dosya: " << (audio.filename.empty() ? "Bilinmiyor" : audio.filename) << "\n";
    cout << "Ornekleme hizi: " << audio.sampleRate << " Hz\n";
    cout << "Kanal sayisi: " << audio.channels << "\n";
    cout << "Toplam ornek: " << audio.totalSamples << "\n";
    cout << "Sure: " << duration << " saniye\n";
    
    // Ses genliÄŸi istatistikleri hesapla
    float maxAmp = 0.0f, avgAmp = 0.0f;
    for (const auto& sample : audio.samples) {
        float absVal = fabs(sample);
        maxAmp = max(maxAmp, absVal);
        avgAmp += absVal;
    }
    avgAmp /= audio.samples.size();
    
    cout << "Maksimum genlik: " << maxAmp << "\n";
    cout << "Ortalama genlik: " << avgAmp << "\n";
    
    // Kanal bilgisi ekle
    if (audio.channels > 1) {
        cout << "Coklu kanal: " << audio.channels << " kanal\n";
        cout << "Toplam ornek sayisi: " << audio.samples.size() << "\n";
        cout << "Her kanal icin ornek sayisi: " << audio.totalSamples << "\n";
    }
}

// Efektler
void echoEffect() {
    if (!audioLoaded) {
        cout << "Lutfen once ses dosyasi yukleyin veya olusturun.\n";
        return;
    }

    int delayMs;
    float decay;
    cout << "Gecikme miktari (ms): ";
    cin >> delayMs;
    cout << "Azalma orani (0-1): ";
    cin >> decay;
    
    // Ornekleme hizina gore gecikme orneklerini hesapla
    int delaySamples = (delayMs * audio.sampleRate) / 1000;
    
    vector<float> result = audio.samples; // Sonuc icin kopya
    
    for (size_t i = delaySamples; i < result.size(); ++i) {
        result[i] += decay * audio.samples[i - delaySamples];
    }
    
    audio.samples = result;
    cout << "[âœ“] Echo efekti uygulandi.\n";
}

void reverbEffect() {
    if (!audioLoaded) {
        cout << "Lutfen once ses dosyasi yukleyin veya olusturun.\n";
        return;
    }

    float decay;
    int numReflections;
    
    cout << "Azalma orani (0-1): ";
    cin >> decay;
    cout << "Yankilama sayisi (1-10): ";
    cin >> numReflections;
    
    numReflections = max(1, min(10, numReflections));
    
    vector<float> result = audio.samples;
    int delaySamples = audio.sampleRate / 10; // 100ms delay
    
    for (int reflection = 1; reflection <= numReflections; reflection++) {
        float currentDecay = pow(decay, reflection);
        int currentDelay = delaySamples * reflection;
        
        for (size_t i = currentDelay; i < result.size(); ++i) {
            result[i] += currentDecay * audio.samples[i - currentDelay];
        }
    }
    
    audio.samples = result;
    cout << "[âœ“] Reverb efekti uygulandi.\n";
}

void distortionEffect() {
    if (!audioLoaded) {
        cout << "Lutfen once ses dosyasi yukleyin veya olusturun.\n";
        return;
    }

    float gain;
    cout << "Distortion miktari (1-20): ";
    cin >> gain;
    
    gain = max(1.0f, min(20.0f, gain));
    
    for (auto& s : audio.samples)
        s = tanh(s * gain);
    cout << "[âœ“] Distortion efekti uygulandi.\n";
}

// Normalleştirme efekti
void normalizeEffect() {
    if (!audioLoaded) {
        cout << "Lutfen once ses dosyasi yukleyin veya olusturun.\n";
        return;
    }
    
    // Maksimum genliÄŸi bul
    float maxAmp = 0.0f;
    for (const auto& sample : audio.samples) {
        float absVal = fabs(sample);
        maxAmp = max(maxAmp, absVal);
    }
    
    // Normalizasyon faktÃ¶rÃ¼
    if (maxAmp > 0.0f) {
        float targetAmp;
        cout << "Hedef maksimum genlik (0-1, genellikle 0.9-0.99): ";
        cin >> targetAmp;
        
        targetAmp = max(0.1f, min(1.0f, targetAmp));
        float factor = targetAmp / maxAmp;
        
        // Tüm Örnekleri normalize et
        for (auto& sample : audio.samples) {
            sample *= factor;
        }
        
        cout << "[âœ“] Normalizasyon efekti uygulandi. Genlik faktoru: " << factor << "\n";
    } else {
        cout << "[-] Normalizasyon yapilamadi - ses cok sessiz.\n";
    }
}

// Fade efekti
void fadeEffect() {
    if (!audioLoaded) {
        cout << "Lutfen once ses dosyasi yukleyin veya olusturun.\n";
        return;
    }
    
    int choice;
    cout << "1) Fade In  2) Fade Out  3) Her Ikisi: ";
    cin >> choice;
    
    int fadeLength;
    cout << "Fade suresi (ms): ";
    cin >> fadeLength;
    
    // ms'yi Ã¶rneklere dÃ¶nÃ¼ÅŸtÃ¼r
    fadeLength = (fadeLength * audio.sampleRate) / 1000;
    fadeLength = min(fadeLength, (int)audio.samples.size() / 4);
    
    if (choice == 1 || choice == 3) {
        // Fade In
        for (int i = 0; i < fadeLength; i++) {
            float factor = (float)i / fadeLength;
            for (int ch = 0; ch < audio.channels; ch++) {
                int idx = i * audio.channels + ch;
                if (idx < audio.samples.size()) {
                    audio.samples[idx] *= factor;
                }
            }
        }
        cout << "[âœ“] Fade In efekti uygulandi.\n";
    }
    
    if (choice == 2 || choice == 3) {
        // Fade Out
        int startPos = audio.samples.size() - fadeLength * audio.channels;
        for (int i = 0; i < fadeLength; i++) {
            float factor = 1.0f - (float)i / fadeLength;
            for (int ch = 0; ch < audio.channels; ch++) {
                int idx = startPos + i * audio.channels + ch;
                if (idx >= 0 && idx < audio.samples.size()) {
                    audio.samples[idx] *= factor;
                }
            }
        }
        cout << "[âœ“] Fade Out efekti uygulandi.\n";
    }
}

// Gercek FFT bazli Spektrum Analizi - YENI!
void spectrumAnalysis() {
    if (!audioLoaded) {
        cout << "Lutfen once ses dosyasi yukleyin veya olusturun.\n";
        return;
    }

    // Analiz iÃ§in ornek boyutunu belirle (2^n olmali)
    const int FRAME_SIZE = 4096; 
    
    // Sadece tek kanalli analize izin ver veya kanallardan birini secin
    vector<float> monoSamples;
    if (audio.channels > 1) {
        cout << "Cok kanalli ses tespit edildi. Hangi kanali analiz etmek istersiniz? (0-" 
             << audio.channels-1 << ", veya -1 tumunu karÄ±ÅŸtÄ±rmak icin): ";
        int channelChoice;
        cin >> channelChoice;
        
        if (channelChoice == -1) {
            // Tum kanallari karistir (downmix)
            monoSamples.resize(audio.totalSamples);
            for (size_t i = 0; i < audio.totalSamples; i++) {
                float mixedSample = 0;
                for (int ch = 0; ch < audio.channels; ch++) {
                    mixedSample += audio.samples[i * audio.channels + ch];
                }
                monoSamples[i] = mixedSample / audio.channels;
            }
        } else {
            // Secilen kanali al
            channelChoice = max(0, min(audio.channels-1, channelChoice));
            monoSamples.resize(audio.totalSamples);
            for (size_t i = 0; i < audio.totalSamples; i++) {
                monoSamples[i] = audio.samples[i * audio.channels + channelChoice];
            }
        }
    } else {
        monoSamples = audio.samples; // Zaten tek kanal
    }
    
    // Ses dosyasinin ortasinda bir frame al
    size_t startIdx = min(monoSamples.size() / 2, monoSamples.size() - FRAME_SIZE);
    
    // FFT iÃ§in kompleks array hazÄ±rla
    vector<complex<double>> fftData(FRAME_SIZE);
    
    // Ornekleri kopyala ve Hanning penceresi uygula
    for (int i = 0; i < FRAME_SIZE; i++) {
        double windowValue = 0.5 * (1 - cos(2 * M_PI * i / (FRAME_SIZE - 1)));
        if (startIdx + i < monoSamples.size()) {
            fftData[i] = monoSamples[startIdx + i] * windowValue;
        } else {
            fftData[i] = 0;
        }
    }
    
    // FFT uygula
    fft(fftData);
    
    // Frekans bantlarini hesapla
    const int NUM_BANDS = 16;
    vector<double> bands(NUM_BANDS);
    vector<float> frequencies(NUM_BANDS);
    
    cout << "\n=== Spektrum Analizi (FFT ile) ===\n";
    
    // Logaritmik frekans bantlarini hesapla
    double minFreq = 20.0; // Hz
    double maxFreq = audio.sampleRate / 2.0; // Nyquist frekansi
    
    for (int band = 0; band < NUM_BANDS; band++) {
        // Logaritmik olarak frekans bantlarini hesapla
        double freqStart = minFreq * pow(maxFreq/minFreq, (double)band / NUM_BANDS);
        double freqEnd = minFreq * pow(maxFreq/minFreq, (double)(band+1) / NUM_BANDS);
        
        // FFT bin'lerini frekans'a dÃ¶nÃ¼ÅŸtÃ¼r
        int binStart = ceil(freqStart * FRAME_SIZE / audio.sampleRate);
        int binEnd = floor(freqEnd * FRAME_SIZE / audio.sampleRate);
        
        binStart = max(1, binStart); // DC bileÅŸenini atla
        binEnd = min(binEnd, FRAME_SIZE/2); // Nyquist'e kadar
        
        // Bu banttaki enerjiyi topla
        double energy = 0;
        for (int bin = binStart; bin <= binEnd; bin++) {
            energy += norm(fftData[bin]);
        }
        
        if (binEnd >= binStart) {
            energy /= (binEnd - binStart + 1); // Ortalama enerji
        }
        
        bands[band] = energy;
        frequencies[band] = (freqStart + freqEnd) / 2; // Orta frekans
    }
    
    // Maksimum deÄŸeri bul ve normalize et
    double maxVal = *max_element(bands.begin(), bands.end());
    if (maxVal > 0) {
        for (auto& val : bands) val = (val / maxVal) * 40; // 40 bar maksimum uzunluk
    }
        
    // Frekans bantlarÄ±nÄ± gÃ¶rselleÅŸtir
    for (int band = 0; band < NUM_BANDS; band++) {
        int numBars = bands[band] + 0.5;
        cout << "[" << (int)frequencies[band] << "Hz]\t";
        for (int j = 0; j < numBars; j++) cout << "#";
        cout << "\n";
    }
}

// Ã‡oklu Kanal KarÄ±ÅŸtÄ±rma - YENI!
void channelMixing() {
    if (!audioLoaded) {
        cout << "Lutfen once ses dosyasi yukleyin veya olusturun.\n";
        return;
    }
    
    if (audio.channels < 2) {
        cout << "[-] Karistirma icin en az 2 kanal gerekli, su anki kanal sayisi: " << audio.channels << "\n";
        return;
    }
    
    cout << "\n=== Coklu Kanal Karistirma ===\n";
    cout << "Mevcut kanal sayisi: " << audio.channels << "\n";
    
    vector<float> channelVolumes(audio.channels, 1.0f);
    vector<int> channelPanning(audio.channels, 0); // -100 (sol) ile 100 (sag) arasi
    
    while (true) {
        cout << "\nKanal Ayarlari:\n";
        for (int i = 0; i < audio.channels; i++) {
            cout << "Kanal " << i+1 << ": Ses = " << channelVolumes[i];
            cout << ", Pan = " << channelPanning[i] << "\n";
        }
        
        cout << "\n1) Kanal ses seviyesi ayarla\n";
        cout << "2) Kanal pan ayarla\n";
        cout << "3) Degisiklikleri uygula ve cik\n";
        cout << "0) Iptal\n";
        cout << "Seciminiz: ";
        
        int choice;
        cin >> choice;
        
        if (choice == 0) return;
        else if (choice == 3) break;
        else if (choice == 1) {
            int channelIdx;
            float volume;
            cout << "Kanal numarasi (1-" << audio.channels << "): ";
            cin >> channelIdx;
            channelIdx--; // 0-bazlÄ± indekse dÃ¶nÃ¼ÅŸtÃ¼r
            
            if (channelIdx < 0 || channelIdx >= audio.channels) {
                cout << "[-] Gecersiz kanal numarasi.\n";
                continue;
            }
            
            cout << "Yeni ses seviyesi (0.0-2.0): ";
            cin >> volume;
            channelVolumes[channelIdx] = max(0.0f, min(2.0f, volume));
        }
        else if (choice == 2) {
            int channelIdx, pan;
            cout << "Kanal numarasi (1-" << audio.channels << "): ";
            cin >> channelIdx;
            channelIdx--; // 0-bazlÄ± indekse dÃ¶nÃ¼ÅŸtÃ¼r
            
            if (channelIdx < 0 || channelIdx >= audio.channels) {
                cout << "[-] Gecersiz kanal numarasi.\n";
                continue;
            }
            
            cout << "Pan degeri (-100=sol, 0=orta, 100=sag): ";
            cin >> pan;
            channelPanning[channelIdx] = max(-100, min(100, pan));
        }
    }
    
    // KarÄ±ÅŸtÄ±rma iÅŸlemini gerÃ§ekleÅŸtir
    // Åu an iÃ§in sadece 2 kanallÄ± stereo destekleyeceÄŸiz (Ã§Ä±ktÄ±)
    int outputChannels = 2; // Stereo Ã§Ä±ktÄ±
    vector<float> mixedAudio(audio.totalSamples * outputChannels);
    
    for (size_t i = 0; i < audio.totalSamples; i++) {
        float leftOut = 0, rightOut = 0;
        
        // Her kanalÄ± kaynak stereo'ya ekle
        for (int ch = 0; ch < audio.channels; ch++) {
            float sample = audio.samples[i * audio.channels + ch] * channelVolumes[ch];
            
            // Pan ayarlamasÄ± yap
            float panValue = channelPanning[ch] / 100.0f; // -1.0 ile 1.0 arasÄ±
            
            // Pan deÄŸeri -1 ise %100 sol, 1 ise %100 saÄŸ, 0 ise %50-%50
            float leftGain = (1.0f - panValue) * 0.5f;
            float rightGain = (1.0f + panValue) * 0.5f;
            
            leftOut += sample * leftGain;
            rightOut += sample * rightGain;
        }
        
        // KarÄ±ÅŸtÄ±rÄ±lmÄ±ÅŸ Ã¶rnekleri Ã§Ä±ktÄ± dizisine yaz
        mixedAudio[i * outputChannels] = leftOut;
        mixedAudio[i * outputChannels + 1] = rightOut;
    }
    
    // Ã‡Ä±ktÄ±yÄ± ana ses verisine kopyala
    audio.samples = mixedAudio;
    audio.channels = outputChannels;
    
    cout << "[âœ“] Coklu kanal karistirma islemi tamamlandi.\n";
    cout << "    Yeni kanal sayisi: " << audio.channels << " (Stereo)\n";
}

// MIDI klavye simulasyonu
void midiSimulator() {
    cout << "\n=== MIDI Klavye Simulasyonu ===\n";
    cout << "Tuslara basin (A: C4, S: D4, D: E4, F: F4, G: G4, H: A4, J: B4)\n";
    cout << "Ust oktavlar icin: W, E, R, T, Y, U\n";
    cout << "ESC ile cikis yapabilirsiniz\n\n";
    
    const int octave = 4; // BaÅŸlangÄ±Ã§ oktavÄ±
    
    // Piyano tuÅŸ frekanslarÄ± (C4'ten baÅŸlayarak)
    vector<float> frequencies = {
        261.63, 293.66, 329.63, 349.23, 392.00, 440.00, 493.88, // C4-B4
        523.25, 587.33, 659.25, 698.46, 783.99, 880.00, 987.77  // C5-B5
    };
    
    while (true) {
        if (_kbhit()) {
            char key = _getch();
            if (key == 27) break; // ESC tuÅŸu
            
            int noteIndex = -1;
            
            // Alt oktav
            switch(tolower(key)) {
                case 'a': noteIndex = 0; break; // C4
                case 's': noteIndex = 1; break; // D4
                case 'd': noteIndex = 2; break; // E4
                case 'f': noteIndex = 3; break; // F4
                case 'g': noteIndex = 4; break; // G4
                case 'h': noteIndex = 5; break; // A4
                case 'j': noteIndex = 6; break; // B4
                
                // Ãœst oktav
                case 'w': noteIndex = 7; break; // C5
                case 'e': noteIndex = 8; break; // D5
                case 'r': noteIndex = 9; break; // E5
                case 't': noteIndex = 10; break; // F5
                case 'y': noteIndex = 11; break; // G5
                case 'u': noteIndex = 12; break; // A5
                case 'i': noteIndex = 13; break; // B5
            }
            
            if (noteIndex >= 0) {
                float freq = frequencies[noteIndex];
                int duration = 300; // ms
                
                // Windows'ta ses Ã§al
                Beep((int)freq, duration);
                
                // Nota bilgisini gÃ¶ster
                string noteName;
                string notes[] = {"C", "D", "E", "F", "G", "A", "B"};
                int noteOctave = (noteIndex < 7) ? octave : octave + 1;
                noteName = notes[noteIndex % 7] + to_string(noteOctave);
                
                cout << "[" << noteName << "] " << freq << " Hz\n";
            }
        }
    }
}

bool isEncryptedUser(const string& username, const string& password) {
    string testUser = decrypt(username);
    string testPass = decrypt(password);

    bool userReadable = all_of(testUser.begin(), testUser.end(), ::isalnum);
    bool passReadable = all_of(testPass.begin(), testPass.end(), ::isalnum);

    return !(userReadable && passReadable);
}



void loadUsersFromFile() {
    users.clear();
    ifstream in("users.txt");

    if (!in || in.peek() == ifstream::traits_type::eof()) {
        cout << "users.txt yok veya bos. Varsayilan kullanicilar ekleniyor...\n";
        users.push_back({ "admin", encrypt("admin123"), true });
        users.push_back({ "user", encrypt("user123"), false });
        saveUsersToFile();
        return;
    }

    string encUser, encPass;
    int isAdmin;

    while (in >> encUser >> encPass >> isAdmin) {
        string username, password;

        if (isEncryptedUser(encUser, encPass)) {
            username = decrypt(encUser);
            password = decrypt(encPass);
        } else {
            username = encUser;
            password = encPass;
        }

        users.push_back({ username, password, static_cast<bool>(isAdmin) });
    }

    in.close();
}



void saveUsersToFile() {
    ofstream out("users.txt");
    for (const auto& user : users) {
        out << encrypt(user.username) << " "
            << encrypt(user.password) << " "
            << user.isAdmin << "\n";
    }
    out.close();
}

// Kullanici yonetimi
void manageUsers() {
    if (!currentUser.isAdmin) {
        cout << "Bu islem sadece admin tarafindan yapilabilir.\n";
        return;
    }

    while (true) {
        cout << "\n=== Kullanici Yonetimi ===\n";
        cout << "1) Tum kullanicilari listele\n";
        cout << "2) Kullanici ekle\n";
        cout << "3) Kullanici sil\n";
        cout << "4) Kullanici yetkisini degistir\n";
        cout << "0) Ana menuye don\n";
        cout << "Seciminiz: ";
        
        int choice;
        cin >> choice;
        
        if (choice == 0) break;
        
        switch (choice) {
            case 1: {
                cout << "\nMevcut Kullanicilar:\n";
                cout << "-------------------\n";
                for (size_t i = 0; i < users.size(); i++) {
                    cout << i+1 << ") " << users[i].username 
                         << " (Yetki: " << (users[i].isAdmin ? "Admin" : "User") << ")\n";
                }
                break;
            }
            case 2: {
                string u, p;
                bool isA;
                cout << "Yeni kullanici adi: "; cin >> u;
                cout << "Sifre: "; cin >> p;
                cout << "Admin mi? (1/0): "; cin >> isA;
                
                // KullanÄ±cÄ± adÄ±nÄ±n benzersiz olup olmadÄ±ÄŸÄ±nÄ± kontrol et
                bool userExists = false;
                for (const auto& user : users) {
                    if (user.username == u) {
                        userExists = true;
                        break;
                    }
                }
                
                if (userExists) {
                    cout << "[-] Bu kullanici adi zaten kullaniliyor!\n";
                } else {
                    users.push_back({u, encrypt(p), isA});
                    saveUsersToFile();

                    cout << "[+] Kullanici eklendi ve kaydedildi.\n";
                }
                break;
            }
            case 3: {
                string u;
                cout << "Silinecek kullanici: "; cin >> u;
                
                // YÃ¶netici kendi hesabÄ±nÄ± silemesin
                if (u == currentUser.username) {
                    cout << "[-] Kendi hesabinizi silemezsiniz!\n";
                    break;
                }
                
                for (auto it = users.begin(); it != users.end(); ++it) {
                    if (it->username == u) {
                        users.erase(it);
                        saveUsersToFile();
                        cout << "[+] Kullanici silindi.\n";
                        break;
                    }
                    
                    if (it == users.end() - 1) {
                        cout << "[-] Kullanici bulunamadi.\n";
                    }
                }
                break;
            }
            case 4: {
                string u;
                cout << "Yetkisi degistirilecek kullanici: "; cin >> u;
                
                for (auto& user : users) {
                    if (user.username == u) {
                        bool newAdminStatus;
                        cout << "Yeni yetki (1: Admin, 0: User): ";
                        cin >> newAdminStatus;
                        
                        user.isAdmin = newAdminStatus;
                        saveUsersToFile();
                        cout << "[+] Kullanici yetkisi guncellendi.\n";
                        break;
                    }
                    
                    if (&user == &users.back()) {
                        cout << "[-] Kullanici bulunamadi.\n";
                    }
                }
                break;
            }
            default:
                cout << "[-] Gecersiz secim.\n";
        }
    }
}

void mainMenu() {
    int secim;
    do {
        cout << "\n=== Dijital Ses Isleme Studyosu ===\n";
        cout << "1) Sahte Ses Uret\n";
        cout << "2) Kaydet (Sifreli)\n";
        cout << "3) Yukle (Sifreli)\n";
        cout << "4) WAV Dosyasi Yukle\n";
        cout << "5) MP3 Dosyasi Yukle \n";
        cout << "6) WAV Dosyasi Kaydet\n";
        cout << "7) Ses Bilgilerini Goster\n";
        cout << "8) Efekt Uygula\n";
        cout << "9) Spektrum Analiz (FFT ile)\n";
        cout << "10) Coklu Kanal Karistirma \n";
        cout << "11) MIDI Klavye\n";
        cout << "12) Kullanici Yonetimi\n";
        cout << "0) Cikis\n";
        cout << "Seciminiz: ";
        cin >> secim;

        switch (secim) {
            case 0: cout << "Cikis yapiliyor...\n"; break;
            case 1: generateDummyAudio(); break;
            case 2: saveAudio(); break;
            case 3: loadAudio(); break;
            case 4: loadWavFile(); break;
            case 5: loadMp3File(); break;
            case 6: saveWavFile(); break;
            case 7: displayAudioInfo(); break;
            case 8: {
                int e;
                cout << "\n=== Efekt Secimi ===\n";
                cout << "1) Echo\n";
                cout << "2) Reverb\n";
                cout << "3) Distortion\n";
                cout << "4) Normalize\n";
                cout << "5) Fade\n";
                cout << "Seciminiz: ";
                cin >> e;
                
                switch (e) {
                    case 1: echoEffect(); break;
                    case 2: reverbEffect(); break;
                    case 3: distortionEffect(); break;
                    case 4: normalizeEffect(); break;
                    case 5: fadeEffect(); break;
                    default: cout << "[-] Gecersiz secim.\n";
                }
                break;
            }
            case 9: spectrumAnalysis(); break;
            case 10: channelMixing(); break;
            case 11: midiSimulator(); break;
            case 12: manageUsers(); break;
            default: cout << "[-] Gecersiz secim.\n";
        }

    } while (secim != 0);
}

int main() {
    srand(time(0));
    SetConsoleTitleA("Dijital Ses Isleme Studyosu v2.0");
    
    // KullanÄ±cÄ±larÄ± yÃ¼kle
    loadUsersFromFile();

    cout << "=== Dijital Ses Studyosu v2.0 Giris ===\n";
    if (!login()) return 0;
    
    mainMenu();
    
    cout << "\nProgrami kapatmak icin bir tusa basiniz...";
    _getch();
    return 0;
}